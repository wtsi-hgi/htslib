#include "htslib/kstring.h"
#include "htslib/bgzf.h"

#include "cram/cram.h"
#include "cram/cram_io.h"
#include "cram/open_trace_file.h"
#include "htslib/hfile.h"

#include <errno.h>
#include <sys/stat.h>

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX
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

/*
 * Populates the ref parameter with a reference to a BGZF file in a cache
 * for the reference genome for a particular M5 string
 *
 * Returns 0 on success
 *        -1 on failure
 */
int m5_to_ref(const char *m5_str, char **name, BGZF *ref) {
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
        expand_cache_path(path, local_cache, m5_str);
        local_path = 1;
    }

#ifndef HAVE_MMAP
    char *path2;
    /* Search local files in REF_PATH; we can open them and return as above */
    if (!local_path && (path2 = find_path(m5_str, ref_path))) {
        strncpy(path, path2, PATH_MAX);
        free(path2);
        if (is_file(path)) // incase it's too long
            local_path = 1;
        // should really error here?
    }
#endif

    if (local_path){
        struct stat sb;

        if(0 == stat(path, &sb) && (ref = bgzf_open(path, "r")))
            /* Found via REF_CACHE or local REF_PATH file */
            return 0;
    }

    /* Look in ref_path */
    if ((mf = open_path_mfile(m5_str, ref_path, NULL))) {
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

        expand_cache_path(path, local_cache, m5_str);
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
            return 0;
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

    int fd = fileno(mf->fp);
    hFILE* hf = hdopen(fd, "r");

    mfdestroy(mf);
    
    ref = bgzf_hopen(hf, "r");
    return 0;
}