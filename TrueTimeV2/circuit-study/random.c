#include "random.h"
#include <stdlib.h>
#include <errno.h>
#ifdef __APPLE__
#include <stdlib.h> // arc4random_buf
#else
#include <sys/random.h>
#endif

int get_true_random_u64(uint64_t *out) {
#ifdef __APPLE__
    arc4random_buf(out, sizeof(*out));
    return 0;
#else
    ssize_t r = getrandom(out, sizeof(*out), 0);
    if (r != sizeof(*out)) return errno;
    return 0;
#endif
}

int get_true_random_double(double *out) {
    uint64_t v;
    int err = get_true_random_u64(&v);
    if (err) return err;
    // top 53 bits => double in [0,1)
    v >>= 11;
    *out = (double)v / (double)(1ULL<<53);
    return 0;
}