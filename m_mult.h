#ifndef MATRIX_MULT_M_MULT_H
#define MATRIX_MULT_M_MULT_H

#include <arm_neon.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

void matrix_multiply(const float *matrix_a, const float *matrix_b, float *matrix_c, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float sum = 0.0f;
            for (int k = 0; k < n; k++) {
                sum += matrix_a[i * n + k] * matrix_b[k * n + j];
            }
            matrix_c[i * n + j] = sum;
        }
    }
}


typedef struct MultiplyArgs MultiplyArgs;
struct MultiplyArgs {
    const float16_t *matrix_a;
    const float16_t *matrix_b;
    float16_t *matrix_c;
    int n;
    int start_row;
    int end_row;
};
void matrix_multiply_simd(const float16_t *matrix_a, const float16_t *matrix_b, float16_t *matrix_c, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(n < 8){
                float sum = 0.0f;
                for (int k = 0; k < n; k++) {
                    sum += matrix_a[i * n + k] * matrix_b[k * n + j];
                }
                matrix_c[i * n + j] = sum;

            }
            else {
                float16x8_t a, b, c = vdupq_n_f16(0);
                for (int k = 0; k < n; k += 8) {
                    a = vld1q_f16(&matrix_a[i * n + k]);
                    for (int l = 0; l < 8; l++) {
                        b[l] = matrix_b[(k + l) * n + j];
                    }
                    c = vfmaq_f16(c, a, b);
                }
                float16x4_t result = vpadd_f16(vget_high_f16(c), vget_low_f16(c));
                result = vpadd_f16(result, vdup_n_f16(0));
                result = vpadd_f16(result, vdup_n_f16(0));
                matrix_c[i * n + j] = vget_lane_f16(result, 0);
            }
        }
    }
}


void matrix_multiply_simd_worker(const MultiplyArgs *args) {
    const float16_t *matrix_a = args->matrix_a;
    const float16_t *matrix_b = args->matrix_b;
    float16_t *matrix_c = args->matrix_c;
    int n = args->n;
    int start_row = args->start_row;
    int end_row = args->end_row;

    for (int i = 0; i < n; i++) {
        for (int j = start_row; j < end_row; j++) {
            if(n < 8){
                float sum = 0.0f;
                for (int k = 0; k < n; k++) {
                    sum += matrix_a[i * n + k] * matrix_b[k * n + j];
                }
                matrix_c[i * n + j] = sum;

            }
            else {
                float16x8_t a, b, c = vdupq_n_f16(0);
                for (int k = 0; k < n; k += 8) {
                    a = vld1q_f16(&matrix_a[i * n + k]);
                    for (int l = 0; l < 8; l++) {
                        b[l] = matrix_b[(k + l) * n + j];
                    }
                    c = vfmaq_f16(c, a, b);
                }
                float16x4_t result = vpadd_f16(vget_high_f16(c), vget_low_f16(c));
                result = vpadd_f16(result, vdup_n_f16(0));
                result = vpadd_f16(result, vdup_n_f16(0));
                matrix_c[i * n + j] = vget_lane_f16(result, 0);
            }
        }
    }
}

void matrix_multiply_simd_threaded(const float16_t *matrix_a, const float16_t *matrix_b, float16_t *matrix_c, int n ) {

    long num_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    if (num_cpus > n) {
        num_cpus = n;
    }


    int rows_per_thread = 1; // Default
    if (num_cpus > 0) {
        rows_per_thread = (int)ceil((double)n / num_cpus);
    }
    if (rows_per_thread == 0 && n > 0) { // If n is small, e.g., n=1, num_cpus=2, ceil(0.5)=1
        rows_per_thread = 1;
    }

    pthread_t threads[num_cpus];
    MultiplyArgs args[num_cpus];
    for (int i = 0; i < num_cpus; i++) {
        args[i].matrix_a = matrix_a;
        args[i].matrix_b = matrix_b;
        args[i].matrix_c = matrix_c;
        args[i].n = n;
        args[i].start_row = i * rows_per_thread;
        args[i].end_row = (i + 1) * rows_per_thread < n ? (i + 1) * rows_per_thread : n;
        pthread_create(&threads[i], NULL, (void *)matrix_multiply_simd_worker, &args[i]);
    }
    for (int i = 0; i < num_cpus; i++) {
        pthread_join(threads[i], NULL);
    }
}
#endif //MATRIX_MULT_M_MULT_H

