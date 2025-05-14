#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define ROW_SIZE 16
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