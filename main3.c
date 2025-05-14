#include <stdio.h>
#include <time.h>
#include "m_mult.h"

#define N 4096
float16_t a[N*N], b[N*N], c[N*N];
int main(){
    struct timespec start, end;
    for(size_t i = 0; i < N * N; i++){
        a[i] = (float) i;
        b[i] = (float) (N * N - i);
    }
    for(size_t i = 1; i <= N; i*=2){
        clock_gettime(CLOCK_MONOTONIC, &start);
        int count = 0;
        RETRY:
        matrix_multiply_simd_threaded(a,b,c,i);
        a[0] = c[0];
        clock_gettime(CLOCK_MONOTONIC, &end);
        long seconds  = end.tv_sec  - start.tv_sec;
        long nanosec  = end.tv_nsec - start.tv_nsec;
        double elapsed = seconds + nanosec*1e-9;
        count++;
        if(elapsed < 0.5){
            goto RETRY;
        }
        elapsed /= count;
        printf("%d,%.17g\n", i, elapsed);
    }
    printf("%f\n", c[0]);
}