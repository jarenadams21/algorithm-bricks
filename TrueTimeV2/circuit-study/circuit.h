#ifndef CIRCUIT_H
#define CIRCUIT_H
#include <complex.h>
#include <stddef.h>

typedef double complex cpx;
struct QuantumRegister { size_t n; cpx *state; };
enum GateType { G_H, G_X, G_CNOT };
struct Gate { enum GateType type; size_t ctrl, targ; };
struct Circuit { size_t n; struct Gate *gates; size_t len, cap; };

int qreg_new(struct QuantumRegister *qr, size_t n);
void qreg_free(struct QuantumRegister *qr);
int qreg_measure_all(struct QuantumRegister *qr, int *bits);

void circ_init(struct Circuit *c, size_t n);
void circ_add_gate(struct Circuit *c, enum GateType t, size_t ctrl, size_t targ);
int circ_run(struct Circuit *c, struct QuantumRegister *qr);
void circ_free(struct Circuit *c);
#endif // CIRCUIT_H