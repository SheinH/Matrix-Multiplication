# Multiplying Matrices

## MIPS Code

The following was the MIPS code I came up with for multiplying two matrices.

```assembly_x86
.data
matrixAFile:    .asciiz "/Users/meow/Downloads/matrixA.bin"
matrixBFile:    .asciiz "/Users/meow/Downloads/matrixB.bin"

bufferSize:     .word   1048576
matrixSize:     .word   128
startTimeLow: .word 0
startTimeHigh: .word 0
endTimeLow:   .word 0
endTimeHigh:  .word 0
bufferA:        .space  1048576
bufferB:        .space  1048576
bufferC:        .space  1048576

.text
                .globl  main



main:
    # load bufferSize into $t2
    lw      $t2,            bufferSize
    lw      $t3,            matrixSize

    # --- Read matrix A ---
    la      $a0,            matrixAFile                     # filename
    la      $a1,            bufferA                         # buffer pointer
    move    $a2,            $t2                             # size
    jal     read_file
    move    $s0,            $v0                             # bytes read from A

    # --- Read matrix B ---
    la      $a0,            matrixBFile
    la      $a1,            bufferB
    move    $a2,            $t2
    jal     read_file
    move    $s1,            $v0                             # bytes read from B

    li      $s0,            0
    li      $s1,            0
    li $v0, 30          # Syscall code for system time
    syscall
    sw $a0, startTimeLow # Store lower 32 bits of start time
    sw $a1, startTimeHigh # Store higher 32 bits of start time
    jal     matrix_mult
    li $v0, 30          # Syscall code for system time
    syscall
    sw $a0, endTimeLow   # Store lower 32 bits of end time
    sw $a1, endTimeHigh  # Store higher 32 bits of end time
    lw $t0, startTimeLow
    lw $t1, endTimeLow
    sub $t0, $t1, $t0   # duration = endTimeLow - startTimeLow
    # exit
    li      $v0,            10
    syscall


read_file:
    move    $t2,            $a1
    move    $t3,            $a2

    # open(filename, O_RDONLY, 0)
    li      $a1,            0
    li      $a2,            0
    li      $v0,            13
    syscall
    move    $t0,            $v0

    # read(fd, buffer, size)
    move    $a0,            $t0
    move    $a1,            $t2
    move    $a2,            $t3
    li      $v0,            14
    syscall
    move    $t1,            $v0

    # close(fd)
    move    $a0,            $t0
    li      $v0,            16
    syscall

    # return bytes read
    move    $v0,            $t1
    jr      $ra

    #s0 = i
    #s1 = j
    #s2 = k
    #s3 = matrix A ptr
    #s4 = matrix B ptr
    #s5 = matrix C ptr
    #s6 = sum
    #s7 = rowSize
matrix_mult:
    addi $sp, $sp, -4       # Make space on stack
    sw   $ra, 0($sp)
    lw $s7, matrixSize #s7 = matrix size
    sll $t9, $s7, 2 #t9 = row size in bytes
    la $s5, bufferC
i_condition:
    beq  $s0, $s7, exit #exit if i = matrix size
    move $s1, $zero #j = 0
j_condition:
    beq $s1, $s7, i_increment #break if j = matrix size
inner_loop_setup:
    mul $t0, $s0, $t7
    la $s3, bufferA($t0) #bufferAPtr = bufferA[rowsize * i]
    mul $t0, $s1, 4
    la $s4, bufferB($t0) #bufferBptr = bufferB[j]
    
    move $s6, $zero
    move $s2, $zero
k_condition:
    beq $s2, $s7, j_increment
inner_loop:
    lw $t0, ($s3)
    lw $t1, ($s4)
    mul $t2, $t0, $t1 #multiply matrix a * matrix b
    add $s6, $s6, $t2 # add to running sum
    addi $s2, $s2, 1 #k += 1
    addi $s3, $s3, 4
    add $s4, $s4, $t9
    j k_condition
i_increment:
    addi $s0, $s0, 1
    j i_condition
j_increment:
    sw $s6, ($s5)
    addi $s5, $s5, 4
    addi $s1, $s1, 1
    j j_condition


exit:
    lw   $ra, 0($sp)
    addi $sp, $sp, 4
    jr   $ra
```

I based this mips code on the following C code:

```c
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
```

In order to simplify the programming process and to avoid having to load in words, I stored all loop variables in the registers.

| C Variable in `matrix_multiply` | MIPS Register in `matrix_mult` | Description                                                                                                   |
|---------------------------------|--------------------------------|---------------------------------------------------------------------------------------------------------------|
| `n` (matrix dimension)          | `$s7`                          | Loaded with `matrixSize`. Referred to as `rowSize` in a comment. Used as the loop boundary.                   |
| `i` (outer loop counter)        | `$s0`                          | Commented as `#s0 = i`. Used to iterate through rows of `matrix_c` and `matrix_a`.                            |
| `j` (middle loop counter)       | `$s1`                          | Commented as `#s1 = j`. Used to iterate through columns of `matrix_c` and rows/columns of `matrix_b`.         |
| `k` (inner loop counter)        | `$s2`                          | Commented as `#s2 = k`. Used to iterate for calculating the dot product for each element of `matrix_c`.       |
| `sum` (float accumulator)       | `$s6`                          | Commented as `#s6 = sum`. Used to accumulate the products `matrix_a[i*n+k] * matrix_b[k*n+j]`.                |
| `matrix_a` (pointer to elements) | `$s3`                         | Commented as `#s3 = matrix A ptr`. Holds the memory address of the current element being accessed from `bufferA` (representing `matrix_a`). It's updated to point to `matrix_a[i*n+k]` within the loops. |
| `matrix_b` (pointer to elements) | `$s4`                         | Commented as `#s4 = matrix B ptr`. Holds the memory address of the current element being accessed from `bufferB` (representing `matrix_b`). It's updated to point to `matrix_b[k*n+j]` within the loops. |
| `matrix_c` (pointer to elements) | `$s5`                         | Commented as `#s5 = matrix C ptr`. Holds the memory address of the current element in `bufferC` (representing `matrix_c`) where the computed `sum` is to be stored. It's updated to point to `matrix_c[i*n+j]`. |


I also made use of system calls in order to load in sample matrix data. In order to actually generate the matrix data, I wrote a separate C program, which simply generated two matrices.
One matrix was an identity matrix that just doubled every element, and the other matrix was every number from 1 to n^2.


### Matrix Generator

```c
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define ROW_SIZE 512

int main() {
    FILE *file = fopen("matrixA.bin", "wb");

    for (uint32_t i = 0; i < ROW_SIZE * ROW_SIZE; ++i) {
        int32_t value = (int32_t) i;
        fwrite(&value, sizeof(int32_t), 1, file);
    }

    fclose(file);
    file = fopen("matrixB.bin", "wb");

    int buf[ROW_SIZE] = {0};
    for (uint32_t i = 0; i < ROW_SIZE; ++i) {
        buf[i] = 2;
        fwrite(buf, sizeof(int), ROW_SIZE, file);
        buf[i] = 0;
    }

    fclose(file);
    return EXIT_SUCCESS;
}
```

To load in the test data, I wrote function `read_file` which takes in a filename, a write address, and a size integer and reads the file to the buffer using syscalls.

Then, I ran my matrix multiplication code. In order to time it, I used another system call to get the current time before and after running matrix multiplication.

| `syscall` Code | C equivalent |
|-|-|
| 13 | `open()` |
| 14 | `read()` |
| 16 | `close()` |
| 30 | `clock_gettime()` | 

By adjusting my `matrixSize` I was able to multiply matrices of different sizes. I was able to multiply up to a 128 by 128 matrix of 32-bit integers before the computation time was too long.

## C Code

### Naive Implementation

Next I wrote my C code for matrix multiplication.

```
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
```

I compiled this code with both `-O0` for no optimizations and `-Ofast -march=native -ffast-math` for maximum optimization.

The optimizations gave roughly double the speed. One change I noticed was the use of `fmadd` in the compiler-optimized version.

### SIMD Optimization

In order to make this even faster, I decided to use SIMD. In order to process as many elements as possible with one instruction, I changed the code to use half-precision floating-point numbers (`float16_t`).

Now, I can load 8 floats from matrix A and b into one SIMD register represented as a `float16x8_t`. Then the function `vfmaq_f16` performs a multiply-add on 8 elements in a single instruction.

After calculating sums for one row, I use `vpadd_f16` to sum up the 8 partial sums. The result is stored in matrix C.

**SIMD code:**

```c
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
```


This led to a speed increase of 5 to 20 times depending on the matrix size. Interestingly, larger matrices had less of a performance gain and I suspect this is due to thermal limits since I'm running this on a laptop.

### Threading

My final optimization was the one that increased the performance the most: making use of multithreading. 


```c
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
```

Here, inside of function `matrix_multiply_simd_threaded`, I used `sysconf(_SC_NPROCESSORS_ONLN)` to determine the optimal of threads. Usually, this works out to be equal to the number of CPU cores the machine has.

I then divided the rows of the matrices and gave each thread a portion of the computation. Because the results of one row do not affect the next, there was no inter-thread synchronization needed.

After spawning all threads, I simply join each of them and then my matrix will be computed. This led to a 30x increase over the unoptimized version.
