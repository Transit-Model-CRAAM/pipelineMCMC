/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-2026 The ECLIPSE authors (CRAAM, Universidade Presbiteriana Mackenzie)
 */
/**
 * @file multithreadCheck.c
 * @brief Simple program to verify OpenMP functionality
 * 
 * This program checks if OpenMP is working correctly by running
 * a parallel section and reporting the number of threads used.
 */

#include <stdio.h>
#include <omp.h>

int main() {
    int nthreads, tid;
    
    /* Fork a team of threads */
    #pragma omp parallel private(tid)
    {
        /* Obtain thread number */
        tid = omp_get_thread_num();
        
        /* Only master thread does this */
        if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }
        
        /* All threads print their thread ID */
        printf("Thread %d is running\n", tid);
    }
    
    /* All threads join master thread and terminate */
    printf("OpenMP is working correctly!\n");
    
    return 0;
}