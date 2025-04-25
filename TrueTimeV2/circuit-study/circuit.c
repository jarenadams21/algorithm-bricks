#include "circuit.h"
#include "random.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Apply Hadamard on qubit q
static void apply_h(struct QuantumRegister *qr, size_t q) {
    size_t dim = 1ULL << qr->n;
    cpx *tmp = malloc(dim * sizeof(cpx));
    memcpy(tmp, qr->state, dim*sizeof(cpx));
    double s = 1.0 / sqrt(2.0);
    for (size_t i = 0; i < dim; ++i) {
        size_t bit = (i >> q) & 1;
        size_t j = i ^ (1ULL << q);
        if (!bit) {
            qr->state[i] = s * (tmp[i] + tmp[j]);
            qr->state[j] = s * (tmp[i] - tmp[j]);
        }
    }
    free(tmp);
}

// Apply Pauli-X on qubit q
static void apply_x(struct QuantumRegister *qr, size_t q) {
    size_t dim = 1ULL << qr->n;
    for (size_t i = 0; i < dim; ++i) {
        if (((i >> q) & 1) == 0) {
            size_t j = i | (1ULL << q);
            cpx tmp = qr->state[i];
            qr->state[i] = qr->state[j];
            qr->state[j] = tmp;
        }
    }
}

// Apply CNOT: ctrl -> targ
static void apply_cnot(struct QuantumRegister *qr, size_t ctrl, size_t targ) {
    size_t dim = 1ULL << qr->n;
    for (size_t i = 0; i < dim; ++i) {
        if (((i >> ctrl) & 1) == 1) {
            size_t j = i ^ (1ULL << targ);
            cpx tmp = qr->state[i];
            qr->state[i] = qr->state[j];
            qr->state[j] = tmp;
        }
    }
}

int qreg_new(struct QuantumRegister *qr, size_t n) {
    size_t dim = 1ULL << n;
    qr->n = n;
    qr->state = calloc(dim, sizeof(cpx));
    if (!qr->state) return -1;
    qr->state[0] = 1.0 + 0.0*I;
    return 0;
}

void qreg_free(struct QuantumRegister *qr) {
    free(qr->state);
    qr->state = NULL;
}

int qreg_measure_all(struct QuantumRegister *qr, int *bits) {
    size_t dim = 1ULL << qr->n;
    double *probs = malloc(dim * sizeof(double));
    if (!probs) return -1;
    double sum = 0;
    for (size_t i = 0; i < dim; ++i) {
        double p = creal(qr->state[i])*creal(qr->state[i]) + cimag(qr->state[i])*cimag(qr->state[i]);
        probs[i] = p;
        sum += p;
    }
    if (sum <= 0) sum = 1;
    for (size_t i = 1; i < dim; ++i) probs[i] += probs[i-1];
    double r;
    if (get_true_random_double(&r)) { free(probs); return -1; }
    r *= sum;
    size_t idx = 0;
    while (idx < dim && probs[idx] < r) idx++;
    free(probs);
    memset(qr->state, 0, dim*sizeof(cpx));
    qr->state[idx] = 1.0 + 0.0*I;
    for (size_t i = 0; i < qr->n; ++i) bits[i] = (idx >> i) & 1;
    return 0;
}

void circ_init(struct Circuit *c, size_t n) {
    c->n = n;
    c->len = 0;
    c->cap = 8;
    c->gates = malloc(c->cap * sizeof(*c->gates));
}

void circ_add_gate(struct Circuit *c, enum GateType t, size_t ctrl, size_t targ) {
    if (c->len == c->cap) {
        c->cap *= 2;
        c->gates = realloc(c->gates, c->cap * sizeof(*c->gates));
    }
    c->gates[c->len++] = (struct Gate){t, ctrl, targ};
}

int circ_run(struct Circuit *c, struct QuantumRegister *qr) {
    for (size_t i = 0; i < c->len; ++i) {
        struct Gate *g = &c->gates[i];
        switch (g->type) {
            case G_H:    apply_h(qr, g->ctrl); break;
            case G_X:    apply_x(qr, g->ctrl); break;
            case G_CNOT: apply_cnot(qr, g->ctrl, g->targ); break;
        }
    }
    return 0;
}

void circ_free(struct Circuit *c) {
    free(c->gates);
    c->gates = NULL;
}