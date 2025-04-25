#ifndef RANDOM_H
#define RANDOM_H
#include <stdint.h>

int get_true_random_u64(uint64_t *out);
int get_true_random_double(double *out);

#endif // RANDOM_H