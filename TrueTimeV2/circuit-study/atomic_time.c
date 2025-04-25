/*
 * Provides:
 *  - true hardware randomness (u64 & double)
 *  - quantum register + basic gates & measurement
 *  - atomic timestamp with entropy jitter
 *  - simple CLI for testing & benchmarking
 */

/* Core Code */
#include "atomic_time.h"
#include "random.h"
#include <time.h>

int quantum_timestamp(uint64_t *out) {
    struct timespec ts;
    if (clock_gettime(CLOCK_REALTIME, &ts)) return -1;
    uint64_t nanos = (uint64_t)ts.tv_sec*1000000000ULL + ts.tv_nsec;
    uint64_t rnd;
    if (get_true_random_u64(&rnd)) return -1;
    *out = nanos ^ (rnd & 0xFFFF);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "random.h"
#include "circuit.h"
#include "atomic_time.h"

static void usage(void){
    fprintf(stderr,
        "Usage:\n"
        "  quantum_sdk_c --mode measure [--qubits N] [--gates LIST]\n"
        "  quantum_sdk_c --mode timestamp\n"
        "  quantum_sdk_c --mode circuit --qubits N --gates LIST\n"
        "  quantum_sdk_c --mode bench [--qubits N]\n"
    );
}

// Parse gates string like "H0,X1,CNOT0:2" into circuit
static void parse_gates(struct Circuit *c, const char *s) {
    char *tok, *str = strdup(s);
    tok = strtok(str, ",");
    while (tok) {
        if (tok[0]=='H') {
            size_t q = atoi(tok+1);
            circ_add_gate(c, G_H, q, 0);
        } else if (tok[0]=='X') {
            size_t q = atoi(tok+1);
            circ_add_gate(c, G_X, q, 0);
        } else if (strncmp(tok,"CNOT",4)==0) {
            char *p = tok+4;
            size_t ctrl = atoi(p);
            while (*p && (*p==':'||*p==',')) ++p;
            size_t targ = atoi(p);
            circ_add_gate(c, G_CNOT, ctrl, targ);
        }
        tok = strtok(NULL, ",");
    }
    free(str);
}

int bench_measure(size_t qubits) {
    struct QuantumRegister qr;
    if (qreg_new(&qr, qubits)) return 1;
    int *bits = malloc(qubits * sizeof(int));
    const int iters = 100000;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < iters; i++) qreg_measure_all(&qr, bits);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double dt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9;
    printf("measure_all( q=%zu ) x%d took %.6fs (%.2fus/op)\n", qubits, iters, dt, dt/iters*1e6);
    free(bits);
    qreg_free(&qr);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { usage(); return 1; }
    const char *mode = argv[2];
    size_t qubits = 3;
    char *gates = NULL;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--qubits") && i+1 < argc) qubits = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--gates") && i+1 < argc) gates = argv[++i];
    }
    if (!strcmp(mode, "measure")) {
        struct Circuit c; circ_init(&c, qubits);
        if (gates) parse_gates(&c, gates);
        struct QuantumRegister qr; qreg_new(&qr, qubits);
        circ_run(&c, &qr);
        int *bits = malloc(qubits * sizeof(int));
        qreg_measure_all(&qr, bits);
        printf("measured:");
        for (size_t i = 0; i < qubits; i++) printf(" %d", bits[i]);
        printf("\n");
        free(bits);
        circ_free(&c); qreg_free(&qr);
    } else if (!strcmp(mode, "timestamp")) {
        uint64_t ts;
        quantum_timestamp(&ts);
        printf("quantum timestamp: %llu\n", (unsigned long long)ts);
    } else if (!strcmp(mode, "circuit")) {
        struct Circuit c; circ_init(&c, qubits);
        if (!gates) { usage(); return 1; }
        parse_gates(&c, gates);
        struct QuantumRegister qr; qreg_new(&qr, qubits);
        circ_run(&c, &qr);
        int *bits = malloc(qubits * sizeof(int));
        qreg_measure_all(&qr, bits);
        printf("circuit measured:");
        for (size_t i = 0; i < qubits; i++) printf(" %d", bits[i]);
        printf("\n");
        free(bits);
        circ_free(&c); qreg_free(&qr);
    } else if (!strcmp(mode, "bench")) {
        return bench_measure(qubits);
    } else {
        usage();
        return 1;
    }
    return 0;
}