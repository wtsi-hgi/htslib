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
#include <string.h>

int main(int argc, char **argv) {
    const char* m5_str = "ac37ec46683600f808cdd41eac1d55cd";
    BGZF* bgzf;
    char* file_path;
    int64_t file_size;

    int error_code = EXIT_SUCCESS;

    if (m5_to_ref(m5_str, &bgzf, &file_size, &file_path) != 0){
        printf("Error in m5_to_ref");
        return EXIT_FAILURE;
    }

    if (strlen(file_path) == 0){
        printf("File path is empty");
        error_code = EXIT_FAILURE;
    }

    if(file_size <= 0){
        printf("Invalid file size '%lli'", file_size);
        error_code = EXIT_FAILURE;
    }

    int code;
    if((code = bgzf_check_EOF(bgzf))){
        printf("Invalid EOF of bgzf. Error code %i", code);
        error_code = EXIT_FAILURE;
    }

    char* bgzf_data = malloc(100);
    ssize_t new_file_size = bgzf_read(bgzf, bgzf_data, 100);

    printf("file sizes: %lli %lli", (long long)new_file_size, (long long)file_size);

    if (strlen(bgzf_data) == 0){
        printf("Invalid file length");
        error_code = EXIT_FAILURE;
    }

    return error_code;
}