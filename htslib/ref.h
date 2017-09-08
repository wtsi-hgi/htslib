/// @file ref.h
/// Reference genome fetching
/*  
    Copyright (C) 2017 Genome Research Ltd.

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
DEALINGS IN THE SOFTWARE.  */


#ifndef HTSLIB_VCF_H
#define HTSLIB_VCF_H

#include "bgzf.h"
#include "htslib/hfile.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    hFILE base;
    hFILE* innerhf;
    int length;
    char* file_name;
    void* mf;
    int is_mem_hfile;
} hFILE_ref;

/* m5_to_ref() - populates the ref parameter with the reference genome 
 * for a given m5 string
 * @param m5_str: The m5 string to query
 * @param ref: This function modifies this parameter with a hFile of the reference
 * genome corresponding to the given m5 string
 * @return 0 on success
 *        -1 on failure
 * 
 * Note: This function is not currently thread safe, so the caller need to
 * acquire locks before calling this
 */
int m5_to_ref(const char *m5_str, hFILE** ref);

#ifdef __cplusplus
}
#endif

#endif
