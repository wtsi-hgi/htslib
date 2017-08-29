/* test/test-ref.c -- ref unit tests

   Copyright (C) 2017 Genome Research Ltd

   Author: Thomas Hickman <th10@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/

#include "htslib/ref.h"
#include "htslib/bgzf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char * mystrncpy(char *dest, const char *src, size_t n)
{
    size_t i;

    for (i = 0; i < n && src[i] != '\0'; i++)
        dest[i] = src[i];
    for ( ; i < n; i++)
        dest[i] = '\0';

    return dest;
}

int main(int argc, char **argv) {
    const char* m5_str = "bbf4de6d8497a119dda6e074521643dc";
    
    int error_code = EXIT_SUCCESS;
    
    Ref ref;
    
    char template[] = "/tmp/htslib_testXXXXXX";
    char* tmp_dir = mkdtemp(template);
    
    if (tmp_dir == NULL){
        printf("Error creating tmp dir\n");
        goto end;
    }
    
    char* tmp = getenv("REF_CACHE");
    char* prev_REF_CACHE;

    if(tmp != NULL){
        int len = strlen(tmp);
        prev_REF_CACHE = malloc(len + 1);
        strcpy(prev_REF_CACHE, tmp);
    }
    else{
        prev_REF_CACHE = NULL;
    }

    setenv("REF_CACHE", template, 1);    
    
    int i;

    for(i = 0;i < 2;i++){
        if (m5_to_ref(m5_str, &ref) != 0){
            printf("Error in m5_to_ref\n");
            error_code = EXIT_FAILURE;
            break;
        }
        
        if(ref.sz <= 0){
            printf("Invalid file size '%lli'\n", (long long)ref.sz);
            error_code = EXIT_FAILURE;
            break;
        }
        
        if (i == 0){
            // Should read from the network
            if(!ref.seq){
                printf("m5_to_ref doesn't populate seq when loading over the network\n");
            }            
        }
        else{
            // Read from the cache
            if(!ref.bgzf){
                printf("When using the cache, m5_to_ref doesn't populate bgzf\n");
                error_code = EXIT_FAILURE;
                goto close_ref;
            }

            if (!ref.name || strlen(ref.name) == 0){
                printf("File path is empty\n");
                error_code = EXIT_FAILURE;
            }
        
            char* bgzf_data = malloc(100);
            ssize_t new_file_size = bgzf_read(ref.bgzf, bgzf_data, 100);
        
            if (new_file_size < 0){
                printf("Invalid file length\n");
                error_code = EXIT_FAILURE;
                goto close_ref;
            }

            free(bgzf_data);
        }

        close_ref:

        if(ref_close(&ref) != 0){
            printf("Failed to close ref\n");
            error_code = EXIT_FAILURE;
        }
    }

    end:
    if(prev_REF_CACHE){
        setenv("REF_CACHE", prev_REF_CACHE, 1);
    }
    
    return error_code;
}