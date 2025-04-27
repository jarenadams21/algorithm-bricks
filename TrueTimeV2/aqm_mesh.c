#define _POSIX_C_SOURCE 200809L
#if defined(__APPLE__)
#include <stdlib.h>  // arc4random_buf
#include <sys/types.h>
extern void arc4random_buf(void *buf, size_t nbytes);
static ssize_t getrandom(void *buf, size_t buflen, unsigned flags) {
    (void)flags;
    arc4random_buf(buf, buflen);
    return (ssize_t)buflen;
}
#else
#include <fcntl.h>
#include <unistd.h>
static ssize_t getrandom(void *buf, size_t buflen, unsigned flags) {
    // Fallback for non-Apple: read from /dev/urandom
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) { perror("open /dev/urandom"); exit(1); }
    ssize_t r = read(fd, buf, buflen);
    close(fd);
    return r;
}
#endif

#include <stdint.h>
#include <inttypes.h>   // for PRIu64
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>

/**
 * Quantum Mesh Vortex Simulation with Spectral-Aware Topology
 * – Enhanced with stakeholder-ready collapse at 40%,
 *   superfluid transition at 90%, and comprehensive phase-space logging.
 */

#define NODES            50
#define MODES            9
#define ROUNDS           30
#define BASE_DT          -1e42
#define STEPS            ROUNDS
#define DIM              5
#define REWIRE_INTERVAL  1
#define BOUNDARY         10.0
#define PITER            30    // power-iteration iterations

// Protein types
//   0 = normal, 1 = beneficial, 2 = invasive (“bad”)
typedef struct {
    double coords[DIM];
    int type;
} Protein;
static Protein Proteins[NODES];

// Data structures
typedef struct { int mode; double ts; } Proposal;
typedef struct { int mode, count; double ts, prev, hash; } Record;
typedef struct { int src, dst; Proposal p; } QMessage;

// Physical constants
static const double h_const  = 6.62607015e-34;
static const double kB_const = 1.380649e-23;
static const double T_const  = 300.0;

// Simulation state
static Proposal *Pending[NODES];
static int       pend_cnt[NODES];
static Record   *Ledger[NODES];
static int       ledg_cnt[NODES];
static double    NodePos[NODES][DIM], NodeVel[NODES][DIM], NodeForce[NODES][DIM];
static double    alpha[NODES], beta[NODES], DT_arr[NODES];
static int       group_type[NODES], mod_exp[NODES];

// Topology & messaging
static int       adj[NODES][NODES];
static uint64_t  msgs_sent[NODES], msgs_recv[NODES];

#define QSIZE 1024
static QMessage qbuf[QSIZE];
static int     qhead = 0, qtail = 0;

// Liveness tracking
static int lastLedgerCnt[NODES] = {0};

// Spectral values
static double spectral_radius = 0.0;
static double spectral_gap    = 0.0;

// Prototypes
static double   hardware_entropy(void);
static double   fDyson(double d);
static double   wCond(int cnt);
static uint64_t modular_pow(uint64_t base, int exp, uint64_t mod);
static void     broadcast(int src, Proposal p);
static void     deliver(void);
static void     consensus_step(int node, int step);
static void     compute_forces(void);
static void     integrate_step(void);
static void     apply_boundaries(int node);
static void     rewire_topology(int step);
static void     compute_spectral(void);
static void     log_spectral(int step);
static void     log_metrics(int step);
static void     log_tape(int step,int node,int prev,int new,double ent,double ts);
static void     log_sensors(int step);
static void     log_msg_stats(int step);
static void     log_extras(int step);
static void     log_invasive(int step);
static void     log_positions(int step);
static void     log_phase_space(int step);
static void     log_hardware(int step);

// ----------------------------------------------------------------------------
// hardware_entropy: [0,1)
static double hardware_entropy(void) {
    uint64_t rnd;
    if (getrandom(&rnd, sizeof rnd, 0) != (ssize_t)sizeof rnd) {
        perror("getrandom");
        exit(1);
    }
    return (double)rnd / (double)UINT64_MAX;
}

// fDyson attenuation
static double fDyson(double d) {
    return exp(-h_const/(kB_const*T_const*(d+1e-12)));
}

// Bose-Einstein occupancy factor
static double wCond(int cnt) {
    double x = cnt % ROUNDS + 1;
    return 1.0/(exp(h_const*x/(kB_const*T_const)) - 1.0);
}

// fast modular exponentiation
static uint64_t modular_pow(uint64_t base, int exp, uint64_t mod) {
    uint64_t res = 1;
    base %= mod;
    while(exp>0) {
        if (exp & 1) res = (res*base) % mod;
        base = (base*base) % mod;
        exp >>= 1;
    }
    return res;
}

// ----------------------------------------------------------------------------
// Hardware test hooks
static double read_cpu_temperature(void) {
    FILE *f = fopen("/sys/class/thermal/thermal_zone0/temp", "r");
    if (!f) return NAN;
    int millideg;
    if (fscanf(f, "%d", &millideg) != 1) { fclose(f); return NAN; }
    fclose(f);
    return millideg / 1000.0;
}

static double measure_rng_latency(void) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    hardware_entropy();
    clock_gettime(CLOCK_MONOTONIC, &t1);
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
}

static double measure_clock_jitter(int samples) {
    double max_diff = 0, diff;
    struct timespec prev, curr;
    clock_gettime(CLOCK_MONOTONIC, &prev);
    for (int i = 0; i < samples; i++) {
        clock_gettime(CLOCK_MONOTONIC, &curr);
        diff = fabs((curr.tv_sec - prev.tv_sec) + (curr.tv_nsec - prev.tv_nsec) * 1e-9);
        if (diff > max_diff) max_diff = diff;
        prev = curr;
    }
    return max_diff;
}

// ----------------------------------------------------------------------------
// log_hardware — CPU temp, RNG latency, clock jitter
static void log_hardware(int step) {
    static bool header = true;
    FILE *f = fopen("hardware_stats.csv", "a");
    if (!f) return;
    if (header) {
        fprintf(f, "step,cpu_temp_c,rng_latency_s,clock_jitter_s\n");
        header = false;
    }
    double temp   = read_cpu_temperature();
    double rng    = measure_rng_latency();
    double jitter = measure_clock_jitter(10);
    fprintf(f, "%d,%.3f,%.6f,%.6f\n", step, temp, rng, jitter);
    fclose(f);
}

// log_positions — record NodePos (x,y,z)
static void log_positions(int step) {
    FILE* f = fopen("position_data.csv","a");
    if (step == 1) fprintf(f, "step,node,x,y,z\n");
    for (int i = 0; i < NODES; i++) {
        fprintf(f, "%d,%d,%.6f,%.6f,%.6f\n",
                step, i,
                NodePos[i][0], NodePos[i][1], NodePos[i][2]);
    }
    fclose(f);
}

// log_phase_space — record (x,y) & (vx,vy)
static void log_phase_space(int step) {
    static bool header = true;
    FILE* f = fopen("phase_space.csv","a");
    if (!f) return;
    if (header) {
        fprintf(f, "step,node,x,y,vx,vy\n");
        header = false;
    }
    for (int i = 0; i < NODES; i++) {
        fprintf(f, "%d,%d,%.6e,%.6e,%.6e,%.6e\n",
                step, i,
                NodePos[i][0], NodePos[i][1],
                NodeVel[i][0], NodeVel[i][1]);
    }
    fclose(f);
}

// rewire_topology with collapse & superfluid
static void rewire_topology(int step) {
    if (step % REWIRE_INTERVAL != 0) return;

    // collapse at 40%: remove 50% of edges
    if (step == (int)(0.4 * STEPS)) {
        double prob = 0.5;
        for (int i = 0; i < NODES; i++) {
            for (int j = 0; j < NODES; j++) {
                if (i != j && ((double)rand() / RAND_MAX) < prob) {
                    adj[i][j] = 0;
                }
            }
        }
        return;
    }

    // superfluid at 90%: full mesh + velocity boost
    if (step == (int)(0.9 * STEPS)) {
        for (int i = 0; i < NODES; i++) {
            for (int j = 0; j < NODES; j++) {
                if (i != j) adj[i][j] = 1;
            }
            for (int k = 0; k < DIM; k++) {
                NodeVel[i][k] *= 1.5;
            }
        }
        return;
    }

    // default invasive-biased random rewiring
    for (int i = 0; i < NODES; i++) {
        for (int j = 0; j < NODES; j++) {
            double p = 0.4;
            if ((Proteins[i].type == 2) ^ (Proteins[j].type == 2)) {
                p = 0.8;
            }
            adj[i][j] = (i != j) && ((double)rand() / RAND_MAX < p);
        }
    }
}

// compute λ₁ and λ₂ of adj via power-iteration + deflation
static void compute_spectral(void) {
    static double x[NODES], y[NODES], v1[NODES], v2[NODES];

    // principal eigenvector
    for (int i = 0; i < NODES; i++) x[i] = 1.0;
    for (int it = 0; it < PITER; it++) {
        for (int i = 0; i < NODES; i++) {
            double sum = 0;
            for (int j = 0; j < NODES; j++) sum += adj[i][j] * x[j];
            y[i] = sum;
        }
        double norm = 0;
        for (int i = 0; i < NODES; i++) norm += y[i] * y[i];
        norm = sqrt(norm);
        for (int i = 0; i < NODES; i++) x[i] = y[i] / norm;
    }
    double lambda1 = 0;
    for (int i = 0; i < NODES; i++) {
        for (int j = 0; j < NODES; j++) {
            lambda1 += x[i] * adj[i][j] * x[j];
        }
        v1[i] = x[i];
    }

    // second eigenvector deflated
    for (int i = 0; i < NODES; i++) v2[i] = (i == 0 ? 0.0 : 1.0);
    double dot = 0;
    for (int i = 0; i < NODES; i++) dot += v1[i] * v2[i];
    for (int i = 0; i < NODES; i++) v2[i] -= dot * v1[i];
    double norm2 = 0;
    for (int i = 0; i < NODES; i++) norm2 += v2[i] * v2[i];
    norm2 = sqrt(norm2);
    for (int i = 0; i < NODES; i++) v2[i] /= norm2;

    for (int it = 0; it < PITER; it++) {
        for (int i = 0; i < NODES; i++) {
            double sum = 0;
            for (int j = 0; j < NODES; j++) sum += adj[i][j] * v2[j];
            y[i] = sum;
        }
        double proj = 0;
        for (int i = 0; i < NODES; i++) proj += y[i] * v1[i];
        for (int i = 0; i < NODES; i++) y[i] -= lambda1 * proj * v1[i];
        norm2 = 0;
        for (int i = 0; i < NODES; i++) norm2 += y[i] * y[i];
        norm2 = sqrt(norm2);
        for (int i = 0; i < NODES; i++) v2[i] = y[i] / norm2;
    }
    double lambda2 = 0;
    for (int i = 0; i < NODES; i++) {
        for (int j = 0; j < NODES; j++) {
            lambda2 += v2[i] * adj[i][j] * v2[j];
        }
    }

    spectral_radius = lambda1;
    spectral_gap    = lambda1 - lambda2;
}

// log_spectral — write λ₁, λ₂, gap
static void log_spectral(int step) {
    FILE *f = fopen("spectral.csv","a");
    if (!f) return;
    if (step == 1) fprintf(f,"step,lambda1,lambda2,spectral_gap\n");
    fprintf(f,"%d,%.6f,%.6f,%.6f\n",
            step,
            spectral_radius,
            spectral_radius - spectral_gap,
            spectral_gap);
    fclose(f);
}

// compute_forces & integrate_step & apply_boundaries
static void compute_forces(void) {
    for (int i = 0; i < NODES; i++)
        for (int k = 0; k < DIM; k++)
            NodeForce[i][k] = 0.0;

    double spec_factor = 1.0 + spectral_gap;
    for (int i = 0; i < NODES; i++) {
        for (int j = i + 1; j < NODES; j++) {
            double dx[DIM], dist2 = 0;
            for (int k = 0; k < DIM; k++) {
                dx[k] = NodePos[j][k] - NodePos[i][k];
                dist2 += dx[k] * dx[k];
            }
            double d = sqrt(dist2) + 1e-12;
            double rep = 1.0/(d*d), att = d*d;
            double w = spec_factor * ((group_type[i] == group_type[j])?1.0:1.5);
            if (Proteins[i].type == 2 || Proteins[j].type == 2) w *= 2.0;
            double net = w * (att - rep) * fDyson(d) * wCond(ledg_cnt[i] + ledg_cnt[j]);
            for (int k = 0; k < 2; k++) {
                double f = net * (dx[k]/d);
                NodeForce[i][k] += f;
                NodeForce[j][k] -= f;
            }
        }
    }
}

static void integrate_step(void) {
    for (int i = 0; i < NODES; i++) {
        for (int k = 0; k < 2; k++) {
            NodeVel[i][k] += DT_arr[i] * NodeForce[i][k];
            NodePos[i][k] += DT_arr[i] * NodeVel[i][k];
        }
        apply_boundaries(i);
    }
}

static void apply_boundaries(int i) {
    for (int k = 0; k < 2; k++) {
        if (NodePos[i][k] > BOUNDARY) {
            NodePos[i][k] = BOUNDARY;
            NodeVel[i][k] *= -1;
        }
        if (NodePos[i][k] < -BOUNDARY) {
            NodePos[i][k] = -BOUNDARY;
            NodeVel[i][k] *= -1;
        }
    }
}

// messaging: broadcast, deliver, consensus_step
static void broadcast(int src, Proposal p) {
    for (int dst = 0; dst < NODES; dst++) {
        if (dst == src || !adj[src][dst]) continue;
        qbuf[qtail % QSIZE] = (QMessage){ src, dst, p };
        qtail++;
        msgs_sent[src]++;
        msgs_recv[dst]++;
    }
}

static void deliver(void) {
    while (qhead < qtail) {
        QMessage m = qbuf[qhead % QSIZE];
        qhead++;
        Pending[m.dst][pend_cnt[m.dst]++] = m.p;
    }
}

static void consensus_step(int node, int step) {
    deliver();
    if (pend_cnt[node] > 0) {
        Proposal pr = Pending[node][0];
        int prev = (ledg_cnt[node] > 0 ? Ledger[node][ledg_cnt[node]-1].mode : -1);
        Record rec = { pr.mode, ledg_cnt[node]+1, pr.ts, (double)prev, 0.0 };
        uint64_t hval = modular_pow(
            (uint64_t)(pr.ts*1e6) ^ (uint64_t)prev,
            mod_exp[node],
            0xFFFFFFFFFFFFULL
        );
        rec.hash = (double)hval;
        Ledger[node][ledg_cnt[node]++] = rec;
        log_tape(step, node, prev, rec.mode, pr.ts, step * DT_arr[node]);
    }
    pend_cnt[node] = 0;
}

// logging: metrics, tape, sensors, msg_stats, extras, invasive
static void log_metrics(int step) {
    int commits = 0;
    for (int i = 0; i < NODES; i++) commits += ledg_cnt[i];
    double sh = 0, li = 0;
    for (int i = 0; i < NODES; i++) {
        double p = 1.0/NODES;
        sh -= p*log(p);
    }
    for (int i = 0; i < NODES; i++) {
        li += NodePos[i][0]*NodeVel[i][1] - NodePos[i][1]*NodeVel[i][0];
    }
    double et = T_const + hardware_entropy()*kB_const;

    FILE* f = fopen("chunk_data.csv","a");
    if (step == 1) fprintf(f,"step,commits,shannon,liouville,eff_temp\n");
    fprintf(f,"%d,%d,%.6e,%.6e,%.6f\n", step, commits, sh, li, et);
    fclose(f);
    log_msg_stats(step);
}

static void log_tape(int step,int node,int prev,int new,double ent,double ts) {
    FILE* f = fopen("tape_data.csv","a");
    if (step==1 && node==0) fprintf(f,"step,node,prev,new,entropy,time\n");
    fprintf(f,"%d,%d,%d,%d,%.6e,%.6e\n", step,node,prev,new,ent,ts);
    fclose(f);
}

static void log_sensors(int step) {
    FILE* f = fopen("sensor_data.csv","a");
    if (step==1) fprintf(f,"step,i,j,heat\n");
    for (int i = 0; i < NODES; i++) {
        for (int j = i+1; j < NODES; j++) {
            double heat = fabs(NodeForce[i][0] - NodeForce[j][0])
                        + fabs(NodeForce[i][1] - NodeForce[j][1]);
            fprintf(f,"%d,%d,%d,%.6e\n", step,i,j,heat);
        }
    }
    fclose(f);
}

static void log_msg_stats(int step) {
    FILE* f = fopen("message_stats.csv","a");
    if (step==1) fprintf(f,"step,node,sent,recv\n");
    for (int i = 0; i < NODES; i++) {
        fprintf(f,"%d,%d,%" PRIu64 ",%" PRIu64 "\n",
                step,i,msgs_sent[i],msgs_recv[i]);
        msgs_sent[i] = msgs_recv[i] = 0;
    }
    fclose(f);
}

static void log_extras(int step) {
    int disagree = 0, totalPairs = NODES*(NODES-1)/2;
    for (int i = 0; i < NODES; i++) {
        int mi = (ledg_cnt[i]>0?Ledger[i][ledg_cnt[i]-1].mode:-1);
        for (int j = i+1; j < NODES; j++) {
            int mj = (ledg_cnt[j]>0?Ledger[j][ledg_cnt[j]-1].mode:-1);
            if (mi != mj) disagree++;
        }
    }
    double safety = (double)disagree/totalPairs;
    int live = 0;
    for (int i = 0; i < NODES; i++) {
        if (ledg_cnt[i] > lastLedgerCnt[i]) {
            live++;
            lastLedgerCnt[i] = ledg_cnt[i];
        }
    }
    double liveness = (double)live/NODES;
    int edges = 0;
    for (int i = 0; i < NODES; i++) {
        for (int j = i+1; j < NODES; j++) {
            if (adj[i][j]) edges++;
        }
    }
    double edge_frac = (double)edges/totalPairs;

    FILE* f = fopen("metrics.csv","a");
    if (step==1) fprintf(f,"step,safety,liveness,edge_frac\n");
    fprintf(f,"%d,%.6f,%.6f,%.6f\n", step,safety,liveness,edge_frac);
    fclose(f);
}

static void log_invasive(int step) {
    FILE* f = fopen("invasive_stats.csv","a");
    if (step==1) fprintf(f,"step,invasive_normal_pairs\n");
    int count = 0;
    for (int i = 0; i < NODES; i++) {
        for (int j = i+1; j < NODES; j++) {
            if (adj[i][j] &&
               ((Proteins[i].type==2 && Proteins[j].type!=2) ||
                (Proteins[j].type==2 && Proteins[i].type!=2)))
            {
                count++;
            }
        }
    }
    fprintf(f,"%d,%d\n", step,count);
    fclose(f);
}

// ----------------------------------------------------------------------------
int main(void) {
    // remove old logs
    const char* files[] = {
        "chunk_data.csv","tape_data.csv","sensor_data.csv",
        "message_stats.csv","metrics.csv","spectral.csv",
        "invasive_stats.csv","position_data.csv",
        "phase_space.csv","hardware_stats.csv"
    };
    for (size_t i = 0; i < sizeof(files)/sizeof(*files); i++)
        remove(files[i]);

    srand((unsigned)time(NULL));

    // Assign protein types (10% invasive)
    for (int i = 0; i < NODES; i++) {
        Proteins[i].type = ((double)rand()/RAND_MAX < 0.1) ? 2 : 0;
    }

    // Initialize state
    for (int i = 0; i < NODES; i++) {
        Pending[i] = calloc(MODES, sizeof *Pending[i]);
        Ledger[i]  = calloc(STEPS, sizeof *Ledger[i]);
        pend_cnt[i] = ledg_cnt[i] = lastLedgerCnt[i] = 0;
        group_type[i] = i % 2;
        mod_exp[i]    = 3 + (i % 5);

        // Initial topology
        for (int j = 0; j < NODES; j++) {
            double p = 0.4;
            if ((Proteins[i].type == 2) ^ (Proteins[j].type == 2)) {
                p = 0.8;
            }
            adj[i][j] = (i != j) && ((double)rand()/RAND_MAX < p);
        }

        alpha[i]  = 0.5 + 0.5*(i/(double)NODES);
        beta[i]   = 0.1 + 0.05*(i/(double)NODES);
        DT_arr[i] = BASE_DT * (1 + 0.2*((double)rand()/RAND_MAX - 0.5));

        for (int k = 0; k < DIM; k++) {
            NodePos[i][k] = 2.0*((double)rand()/RAND_MAX) - 1.0;
            NodeVel[i][k] = 0.0;
        }

        msgs_sent[i] = msgs_recv[i] = 0;
    }

    // Main simulation loop
    for (int step = 1; step <= STEPS; step++) {
        rewire_topology(step);
        compute_spectral();    log_spectral(step);

        for (int i = 0; i < NODES; i++) {
            Proposal p = { rand()%MODES, hardware_entropy() };
            broadcast(i, p);
        }
        for (int i = 0; i < NODES; i++) {
            consensus_step(i, step);
        }

        compute_forces();
        integrate_step();

        log_metrics(step);
        log_sensors(step);
        log_extras(step);
        log_invasive(step);
        log_positions(step);
        log_phase_space(step);
        log_hardware(step);
    }

    return 0;
}