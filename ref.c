/*
Copyright (c) 2012-2017 Genome Research Ltd.
Authors: James Bonfield <jkb@sanger.ac.uk>, Thomas Hickman <th10@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "htslib/kstring.h"
#include "htslib/bgzf.h"

#include "cram/cram.h"
#include "cram/cram_io.h"
#include "cram/os.h"
#include "cram/misc.h"
#include "cram/open_trace_file.h"
#include "htslib/hfile.h"
#include "htslib/ref.h"
#include "hfile_internal.h"

#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>

#include <config.h>

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX // was 1024 in open_trace_file
#endif


/*
 * Return the cache directory to use, based on the first of these
 * environment variables to be set to a non-empty value.
 */
static const char *get_cache_basedir(const char **extra)
{
    char *base;

    *extra = "";

    base = getenv("XDG_CACHE_HOME");
    if (base && *base)
        return base;

    base = getenv("HOME");
    if (base && *base)
    {
        *extra = "/.cache";
        return base;
    }

    base = getenv("TMPDIR");
    if (base && *base)
        return base;

    base = getenv("TEMP");
    if (base && *base)
        return base;

    return "/tmp";
}

/*
 * Converts a directory and a filename into an expanded path, replacing %s
 * in directory with the filename and %[0-9]+s with portions of the filename
 * Any remaining parts of filename are added to the end with /%s.
 */
void expand_cache_path(char *path, char *dir, char *fn)
{
    char *cp;

    while ((cp = strchr(dir, '%')))
    {
        strncpy(path, dir, cp - dir);
        path += cp - dir;

        if (*++cp == 's')
        {
            strcpy(path, fn);
            path += strlen(fn);
            fn += strlen(fn);
            cp++;
        }
        else if (*cp >= '0' && *cp <= '9')
        {
            char *endp;
            long l;

            l = strtol(cp, &endp, 10);
            l = MIN(l, strlen(fn));
            if (*endp == 's')
            {
                strncpy(path, fn, l);
                path += l;
                fn += l;
                *path = 0;
                cp = endp + 1;
            }
            else
            {
                *path++ = '%';
                *path++ = *cp++;
            }
        }
        else
        {
            *path++ = '%';
            *path++ = *cp++;
        }
        dir = cp;
    }
    strcpy(path, dir);
    path += strlen(dir);
    if (*fn && path[-1] != '/')
        *path++ = '/';
    strcpy(path, fn);
}

/*
 * Return an integer representation of pthread_self().
 */
static unsigned get_int_threadid()
{
    pthread_t pt = pthread_self();
    unsigned char *s = (unsigned char *)&pt;
    size_t i;
    unsigned h = 0;
    for (i = 0; i < sizeof(pthread_t); i++)
        h = (h << 5) - h + s[i];
    return h;
}

/*
 * Make the directory containing path and any prefix directories.
 */
void mkdir_prefix(char *path, int mode)
{
    char *cp = strrchr(path, '/');
    if (!cp)
        return;

    *cp = 0;
    if (is_directory(path))
    {
        *cp = '/';
        return;
    }

    if (mkdir(path, mode) == 0)
    {
        chmod(path, mode);
        *cp = '/';
        return;
    }

    mkdir_prefix(path, mode);
    mkdir(path, mode);
    chmod(path, mode);
    *cp = '/';
}

int ref_close(hFILE* hf){
    int error_code = 0;
    hFILE_ref *hfref = (hFILE_ref*)hf;
    free(hfref->file_name);

    if(hfref->mf)
        error_code |= mfclose(hfref->mf);

    error_code |= hfref->innerhf->backend->close(hfref->innerhf);

    return error_code;
}

static ssize_t ref_read(hFILE *fp, void *bufferv, size_t nbytes) {
    hFILE_ref *hfref = (hFILE_ref*)fp;
    return hfref->innerhf->backend->read(hfref->innerhf, bufferv, nbytes);
}

static ssize_t ref_write(hFILE *fp, const void *bufferv, size_t nbytes) {
    hFILE_ref *hfref = (hFILE_ref*)fp;
    return hfref->innerhf->backend->write(hfref->innerhf, bufferv, nbytes);
}

static off_t ref_seek(hFILE *fp, off_t offset, int whence) {
    hFILE_ref *hfref = (hFILE_ref*)fp;
    return hfref->innerhf->backend->seek(hfref->innerhf, offset, whence);
}

static int ref_flush(hFILE *fp) {
    hFILE_ref *hfref = (hFILE_ref*)fp;
    if(hfref->innerhf->backend->flush)
        return hfref->innerhf->backend->flush(hfref->innerhf);
    else
        return 0;
}

static const struct hFILE_backend ref_backend =
{
    ref_read, ref_write, ref_seek, ref_flush, ref_close
};

/*
 * Modifies the base_hFILE parameter to contain the length and filename information
 * as a hFILE_ref
*/
static hFILE* open_hfile_ref(hFILE* base_hFILE, char* file_name, int length, mFILE* mf, int is_mem_hfile){
    hFILE_ref* hfref = (hFILE_ref *) hfile_init(sizeof (hFILE_ref), "r", 0);

    hfref->file_name = file_name;
    hfref->length = length;
    hfref->mf = mf;
    hfref->innerhf = base_hFILE;
    hfref->base.backend = &ref_backend;
    hfref->is_mem_hfile = is_mem_hfile;

    return &hfref->base;
}

#define READ_CHUNK_SIZE 8192

/*
 * Puts the char* into a hfile
*/
hFILE* mFILE_to_hFILE_ref(mFILE* mf, char* path){
    hFILE *hf = NULL;
    size_t sz;
    char* seq;
    seq = mfsteal(mf, &sz);

    if(!seq){
        seq = mf->data;
    }
    else{
        mf = NULL;
    }
    
    hf = hopen("mem:", "r:", seq, sz);

    return open_hfile_ref(hf, path, sz, mf, 1);
}

// Public functions

int m5_to_ref(const char *m5_str, hFILE** ref) {
    char *ref_path = getenv("REF_PATH");
    char *path = malloc(PATH_MAX); // may be exposed to the outside
    char path_tmp[PATH_MAX];
    char cache[PATH_MAX], cache_root[PATH_MAX];
    char *local_cache = getenv("REF_CACHE");
    hFILE *hf;

    cache_root[0] = '\0';

    int local_path = 0;

    if (!ref_path || *ref_path == '\0') {
        /*
        * If we have no ref path, we use the EBI server.
        * However to avoid spamming it we require a local ref cache too.
        */

        ref_path = "http://www.ebi.ac.uk:80/ena/cram/md5/%s";
        if (!local_cache || *local_cache == '\0') {
            const char *extra;
            const char *base = get_cache_basedir(&extra);
            snprintf(cache_root, PATH_MAX, "%s%s/hts-ref", base, extra);
            snprintf(cache, PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, extra);
            local_cache = cache;
            hts_log_info("Populating local cache: %s", local_cache);
        }
    }

    /* Use cache if available */
    if (local_cache && *local_cache) {
        expand_cache_path(path, local_cache, (char *)m5_str);
        local_path = 1;
    }

#ifndef HAVE_MMAP
    char *path2;
    /* Search local files in REF_PATH; we can open them and return as above */
    if (!local_path && (path2 = find_path((char *)m5_str, ref_path))) {
        strncpy(path, path2, PATH_MAX);
        free(path2);
        if (is_file(path)) // incase it's too long
            local_path = 1;
    }
#endif

    if (local_path){
        struct stat sb;
        hFILE* innerhf;

        if((0 == stat(path, &sb)) && (innerhf = hopen(path, "r"))){
            /* Found via REF_CACHE or local REF_PATH file */
            *ref = open_hfile_ref(innerhf, path, sb.st_size, NULL, 0);
            
            return 0;
        }
    }

    char *resolved_file;
    mFILE* mf;
    /* Look in ref_path */
    if (!(mf = open_path_mfile((char*)m5_str, ref_path, NULL, &resolved_file))) {
        hts_log_info("Failed to fetch file. REF_PATH: '%s', M5: '%s'", ref_path, m5_str);
        return -1;
    }

    *ref = mFILE_to_hFILE_ref(mf, resolved_file);

    // open the hfile separately for testing
    hf = hopen(resolved_file, "r");

    /* Populate the local disk cache if required */
    if (local_cache && *local_cache) {
        int pid = (int)getpid();
        unsigned thrid = get_int_threadid();
        hFILE *fp;

        if (*cache_root && !is_directory(cache_root)) {
            hts_log_warning("Creating reference cache directory %s\n"
                            "This may become large; see the samtools(1) manual page REF_CACHE discussion",
                            cache_root);
        }

        expand_cache_path(path, local_cache, (char *)m5_str);
        hts_log_info("Writing cache file '%s'", path);
        mkdir_prefix(path, 01777);

        do {
            // Attempt to further uniquify the temporary filename
            unsigned t = ((unsigned)time(NULL)) ^ ((unsigned)clock());
            thrid++; // Ensure filename changes even if time/clock haven't

            sprintf(path_tmp, "%s.tmp_%d_%u_%u", path, pid, thrid, t);
            fp = hopen(path_tmp, "wx");
        } while (fp == NULL && errno == EEXIST);
        if (!fp) {
            perror(path_tmp);

            // Not fatal - we have the data already so keep going.
        }

        // Stream the file into the cache and check the md5
        hts_md5_context *md5;
        char unsigned md5_buf1[16];
        char md5_buf2[33];

        if (!(md5 = hts_md5_init())) {
            hclose_abruptly(fp);
            unlink(path_tmp);
            
            return -1;
        }

        int read_length;
        char buf[READ_CHUNK_SIZE];

        while ((read_length = hread(hf, buf, READ_CHUNK_SIZE)) > 0) {
            hts_md5_update(md5, buf, read_length);
            if(hwrite(fp, buf, read_length) != read_length){
                perror(path);
                hclose_abruptly(hf);
                return -1;
            }
        }
        
        hts_md5_final(md5_buf1, md5);
        hts_md5_destroy(md5);
        hts_md5_hex(md5_buf2, md5_buf1);

        if (strncmp(m5_str, md5_buf2, 32) != 0) {
            hts_log_error("Mismatching md5sum for downloaded reference");
            hclose_abruptly(fp);
            unlink(path_tmp);
            return -1;
        }
        if (hclose(fp) < 0)
		{
            unlink(path_tmp);
            
            return -1; // Can't continue if we don't have a file handle
		}
		else
		{
			if (0 == chmod(path_tmp, 0444))
				rename(path_tmp, path);
			else{
                unlink(path_tmp);
                return -1;
            }
        }
    }

    free(path); // not used in exported hfile, so can free

    return 0;
}