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

## Performance Data & Analysis

| Square Matrix Size | MIPS (int32) | Unoptimized C | Compiler Optimized | SIMD     | SIMD + Threads | Compiler Increase | SIMD Increase | Threading Increase |
| ------------------ | ------------ | ------------- | ------------------ | -------- | -------------- | ----------------- | ------------- | ------------------ |
| 1                  | 0            | 0.000000      | 0.000000           | 0.000000 | 0.000017       | 1.1825            | 1.1717        | 0.0018             |
| 2                  | 1            | 0.000000      | 0.000000           | 0.000000 | 0.000033       | 1.7253            | 1.6808        | 0.0013             |
| 4                  | 6            | 0.000000      | 0.000000           | 0.000000 | 0.000046       | 6.4945            | 4.6842        | 0.0051             |
| 8                  | 38           | 0.000002      | 0.000000           | 0.000000 | 0.000104       | 9.9204            | 12.8704       | 0.0163             |
| 16                 | 303          | 0.000014      | 0.000001           | 0.000001 | 0.000125       | 13.3314           | 20.0990       | 0.1116             |
| 32                 | 1854         | 0.000114      | 0.000007           | 0.000005 | 0.000119       | 15.7041           | 23.3563       | 0.9555             |
| 64                 | 13824        | 0.000819      | 0.000074           | 0.000045 | 0.000125       | 11.0722           | 18.0742       | 6.5794             |
| 128                | 109788       | 0.006889      | 0.000753           | 0.000589 | 0.000210       | 9.1534            | 11.7017       | 32.8473            |
| 256                |              | 0.056966      | 0.010521           | 0.005098 | 0.001073       | 5.4144            | 11.1733       | 53.1122            |
| 512                |              | 0.459662      | 0.084347           | 0.084100 | 0.014942       | 5.4496            | 5.4657        | 30.7626            |
| 1024               |              | 3.746641      | 0.872080           | 0.685130 | 0.115436       | 4.2962            | 5.4685        | 32.4564            |
| 2048               |              | 60.1373       | 25.8772            | 22.4906  | 2.0150         | 2.3240            | 2.6739        | 29.8454            |
| 4096               |              |               |                    | 216.5711 | 76.6781        |                   |               |                    |


![Comparison of Matrix Multiplication Methods](https://raw.githubusercontent.com/SheinH/Matrix-Multiplication/refs/heads/main/performance_chart.png)

From this data, a few things stand out:

* Hand written SIMD beats out compiler auto-vectorization, though this may be due to using fp16
* SIMD has less of an edge as matrix size increases: this may be due to thermal constraints on the CPU however
* Threading is actually slower than just SIMD for small matrices. It only becomes better once we are dealing with 64x64 matrices or larger.

### Takeaways

1. Handwritten SIMD can still offer a performance edge compared to compiler-generated auto-vectorized code.
2. Threading can offer massive performance increases, but...
3. Every thread spawned has a performance overhead. Threading should only be employed when data sizes are large enough to warrant them.
4. In comparison, SIMD is faster than non SIMD code almost always.
5. Using half-precision floating points can be advantageous since we can fit more individual data points into a SIMD register. If I used single precision floating points, I would be able to perform only half as many calculations per instruction.
