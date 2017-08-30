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
#include "htslib/hfile.h"
#include "htslib/ref.h"
#include "htslib/hfile.h"

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

// copied from open_trace_file

// Original file from open_trace_file

/*
Author: James Bonfield

Copyright (c) 2000-2001 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF 
MOLECULAR BIOLOGY nor the names of its contributors may be used to endorse or 
promote products derived from this software without specific prior written 
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * Tokenises the search path splitting on colons (unix) or semicolons
 * (windows).
 * We also  explicitly add a "./" to the end of the search path
 *
 * Returns: A new search path with items separated by nul chars. Two nul
 *          chars in a row represent the end of the tokenised path.
 * Returns NULL for a failure.
 *
 * The returned data has been malloced. It is up to the caller to free this
 * memory.
 */
char *tokenise_search_path(const char *searchpath) {
    char *newsearch;
    unsigned int i, j;
    size_t len;
#ifdef _WIN32
    char path_sep = ';';
#else
    char path_sep = ':';
#endif

    if (!searchpath)
	searchpath="";

    newsearch = (char *)malloc((len = strlen(searchpath))+5);
    if (!newsearch)
	return NULL;

    for (i = 0, j = 0; i < len; i++) {
	/* "::" => ":". Used for escaping colons in http://foo */
	if (i < len-1 && searchpath[i] == ':' && searchpath[i+1] == ':') {
	    newsearch[j++] = ':';
	    i++;
	    continue;
	}

	/* Handle http:// and ftp:// too without :: */
	if (path_sep == ':') {
	    if ((i == 0 || (i > 0 && searchpath[i-1] == ':')) &&
		(!strncmp(&searchpath[i], "http:",     5) ||
		 !strncmp(&searchpath[i], "ftp:",      4) ||
		 !strncmp(&searchpath[i], "|http:",    6) ||
		 !strncmp(&searchpath[i], "|ftp:",     5) ||
		 !strncmp(&searchpath[i], "URL=http:", 9) ||
		 !strncmp(&searchpath[i], "URL=ftp:",  8))) {
		do {
		    newsearch[j++] = searchpath[i];
		} while (i<len && searchpath[i++] != ':');
		if (searchpath[i] == ':')
		    i++;
		if (searchpath[i]=='/')
		    newsearch[j++] = searchpath[i++];
		if (searchpath[i]=='/')
		    newsearch[j++] = searchpath[i++];
		// Look for host:port
		do {
		    newsearch[j++] = searchpath[i++];
		} while (i<len && searchpath[i] != ':' && searchpath[i] != '/');
		newsearch[j++] = searchpath[i++];
		if (searchpath[i] == ':')
		    i++;
	    }
	}

	if (searchpath[i] == path_sep) {
	    /* Skip blank path components */
	    if (j && newsearch[j-1] != 0)
		newsearch[j++] = 0;
	} else {
	    newsearch[j++] = searchpath[i];
	}
    }

    if (j)
	newsearch[j++] = 0;
    newsearch[j++] = '.';
    newsearch[j++] = '/';
    newsearch[j++] = 0;
    newsearch[j++] = 0;
    
    return newsearch;
}

/*
 * Takes a dirname possibly including % rules and appends the filename
 * to it.
 *
 * Returns expanded pathname or NULL for malloc failure.
 */
static char *expand_path(const char *file, const char *dirname) {
    size_t len = strlen(dirname);
    size_t lenf = strlen(file);
    char *cp, *path;

    path = malloc(len+lenf+2); // worst expansion DIR/FILE
    if (!path)
	return NULL;

    if (dirname[len-1] == '/')
	len--;

    /* Special case for "./" or absolute filenames */
    if (*file == '/' || (len==1 && *dirname == '.')) {
	sprintf(path, "%s", file);
    } else {
	/* Handle %[0-9]*s expansions, if required */
	char *path_end = path;
	*path = 0;
	while ((cp = strchr(dirname, '%'))) {
	    char *endp;
	    long l = strtol(cp+1, &endp, 10);
	    if (*endp != 's') {
		strncpy(path_end, dirname, (endp+1)-dirname);
		path_end += (endp+1)-dirname;
		dirname = endp+1;
		continue;
	    }
	    
	    strncpy(path_end, dirname, cp-dirname);
	    path_end += cp-dirname;
	    if (l) {
		strncpy(path_end, file, l);
		path_end += MIN(strlen(file), l);
		file     += MIN(strlen(file), l);
	    } else {
		strcpy(path_end, file);
		path_end += strlen(file);
		file     += strlen(file);
	    }
	    len -= (endp+1) - dirname;
	    dirname = endp+1;
	}
	strncpy(path_end, dirname, len);
	path_end += MIN(strlen(dirname), len);
	*path_end = 0;
	if (*file) {
	    *path_end++ = '/';
	    strcpy(path_end, file);
	}
    }

    //fprintf(stderr, "*PATH=\"%s\"\n", path);
    return path;
}

/*
 * Searches for file in the directory 'dirname'. If it finds it, it opens
 * it. This also searches for compressed versions of the file in dirname
 * too.
 *
 * Returns mFILE pointer if found
 *         NULL if not
 */
static hFILE *find_file_dir(const char *file, const char *dirname) {
    char *path;
    hFILE *hf = NULL;

    path = expand_path(file, dirname);

    if (is_file(path))
	hf = hopen(path, "rbm");

    free(path);
    return hf;
}

hFILE *find_file_url(const char *file, const char *url) {
    char buf[8192], *cp;
    mFILE *mf = NULL;
    int maxlen = 8190 - strlen(file), len;
    hFILE *hf;

    // FIXME: should this use expand_path
    /* Expand %s for the trace name */
    for (cp = buf; *url && cp - buf < maxlen; url++) {
	if (*url == '%' && *(url+1) == 's') {
	    url++;
	    cp += strlen(strcpy(cp, file));
	} else {
	    *cp++ = *url;
	}
    }
    *cp++ = 0;

    return hopen(buf, "r");
}

/*
 * ------------------------------------------------------------------------
 * Public functions below.
 */

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in 'path' (which is a colon separated list).
 * If 'path' is NULL it uses the RAWDATA environment variable instead.
 *
 * Returns a hFILE pointer when found.
 *           NULL otherwise.
 */
hFILE *open_path_hfile(const char *file, const char *path) {
    char *newsearch;
    char *ele;
    hFILE *fp;

    /* Use path first */
    if (!path)
	path = getenv("RAWDATA"); // FIXME: is this necessary?
    if (NULL == (newsearch = tokenise_search_path(path)))
	return NULL;
    
    /*
     * Step through the search path testing out each component.
     * We now look through each path element treating some prefixes as
     * special, otherwise we treat the element as a directory.
     */
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
	    char *ele2;

        /*
        * '|' prefixing a path component indicates that we do not
        * wish to perform the compression extension searching in that
        * location.
        *
        * NB: this has been removed from the htslib implementation.
        */
        if (*ele == '|') {
            ele2 = ele+1;
        } else {
            ele2 = ele;
        }

        if (0 == strncmp(ele2, "URL=", 4)) {
            if ((fp = find_file_url(file, ele2+4))) {
            free(newsearch);
            return fp;
            }
        } else if (!strncmp(ele2, "http:", 5) ||
            !strncmp(ele2, "ftp:", 4)) {
            if ((fp = find_file_url(file, ele2))) {
            free(newsearch);
            return fp;
            }
        } else if ((fp = find_file_dir(file, ele2))) {
            free(newsearch);
            return fp;
        } 
    }

    free(newsearch);

    return NULL;
}


/*
 * As per open_path_mfile, but searching only for local filenames.
 * This is useful as we may avoid doing a full mfopen and loading
 * the entire file into memory.
 *
 * Returns the expanded pathname if found.
 *         NULL if not
 */
char *find_path(char *file, char *path) {
    char *newsearch;
    char *ele;
    char *outpath = NULL;

    /* Use path first */
    if (!path)
	path = getenv("RAWDATA");
    if (NULL == (newsearch = tokenise_search_path(path)))
	return NULL;
    
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
	char *ele2 = (*ele == '|') ? ele+1 : ele;

	if (!strncmp(ele2, "URL=", 4) ||
	    !strncmp(ele2, "http:", 5) ||
	    !strncmp(ele2, "ftp:", 4)) {
	    continue;
	} else {
	    outpath = expand_path(file, ele2);
	    if (is_file(outpath)) {
		free(newsearch);
		return outpath;
	    } else {
		free(outpath);
	    }
	} 
    }

    free(newsearch);

    return NULL;
}


// Public functions

int ref_close(Ref* ref){
    int fail = 0;

    if (ref->bgzf){
        return bgzf_close(ref->bgzf);
    }
    else {
        free(ref->seq);

        if(ref->mf == NULL){
            return 0;
        }
        else{
            return mfclose(ref->mf);
        }
    }
}

int m5_to_ref(const char *m5_str, Ref* ref) {
    char *ref_path = getenv("REF_PATH");
    char path[PATH_MAX], path_tmp[PATH_MAX];
    char cache[PATH_MAX], cache_root[PATH_MAX];
    char *local_cache = getenv("REF_CACHE");
    mFILE *mf;

    cache_root[0] = '\0';    

    // If ref_path is not defined or some other situation which I don't think happens
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

        if((0 == stat(path, &sb)) && (ref->bgzf = bgzf_open(path, "r"))){
            /* Found via REF_CACHE or local REF_PATH file */
            
            ref->seq = NULL;
            ref->sz = sb.st_size;
            ref->name = path;
            return 0;
        }
    }

    hFILE* hf = open_path_hfile(m5_str, ref_path);

    /* Look in ref_path */
    if (!(mf = open_path_hfile((char *)m5_str, ref_path))) {
        hts_log_error("Failed to fetch file. REF_PATH: '%s', M5: '%s'", ref_path, m5_str);
        return -1;
    }

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

        // Check md5sum
        hts_md5_context *md5;
        char unsigned md5_buf1[16];
        char md5_buf2[33];

        if (!(md5 = hts_md5_init())) {
            hclose_abruptly(fp);
            unlink(path_tmp);
            return -1;
        }
        hts_md5_update(md5, mf->data, mf->size);
        hts_md5_final(md5_buf1, md5);
        hts_md5_destroy(md5);
        hts_md5_hex(md5_buf2, md5_buf1);

        if (strncmp(m5_str, md5_buf2, 32) != 0) {
            hts_log_error("Mismatching md5sum for downloaded reference");
            hclose_abruptly(fp);
            unlink(path_tmp);
            return -1;
        }

        if (hwrite(fp, mf->data, mf->size) != mf->size) {
            perror(path);
        }
        if (hclose(fp) < 0) {
            unlink(path_tmp);
        }
        else {
            if (0 == chmod(path_tmp, 0444))
                rename(path_tmp, path);
            else
                unlink(path_tmp);
        }
    }

    size_t sz;
    
    ref->seq = mfsteal(mf, &sz);
    if (ref->seq)
    {
        ref->mf = NULL;
    }
    else
    {
        // keep mf around as we couldn't detach
        ref->seq = mf->data;
        ref->mf = mf;
    }
    ref->bgzf = NULL;
    ref->name = NULL;

    return 0;
}