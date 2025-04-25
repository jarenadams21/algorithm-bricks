/* Atomic Time Service */
#ifndef ATOMIC_TIME_H
#define ATOMIC_TIME_H
#include <stdint.h>

int quantum_timestamp(uint64_t *out);

#endif // ATOMIC_TIME_H