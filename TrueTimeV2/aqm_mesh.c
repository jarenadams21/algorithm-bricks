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
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

/**
 * Quantum Mesh Vortex Simulation with per-node and dynamic topology
 * - Nested vortex shapes via per-node α, β radii and Lissajous angles
 * - Randomized DT per node for braided timing
 * - Time-varying adjacency every REWIRE_INTERVAL
 * - Boundary conditions: bounce in square [-1,+1]^2
 * Outputs: chunk_data.csv, tape_data.csv, sensor_data.csv
 */

#define NODES            3
#define MODES            1
#define ROUNDS           1
#define BASE_DT          10e-42
#define STEPS            ROUNDS
#define DIM              5
#define PCA_DIM          3
#define REWIRE_INTERVAL  1
#define BOUNDARY         1.0

// Data structures
typedef struct { int mode; double ts; } Proposal;
typedef struct { int mode, count; double ts, prev, hash; } Record;
typedef struct { int src, dst; Proposal p; } QMessage;

// Physical constants
static const double h_const  = 6.62607015e-34;
static const double kB_const = 1.380649e-23;
static const double T_const  = 30.0;
static const double Lambda   = 1.0;

// Simulation state
static Proposal *Pending[NODES];
static int       pend_cnt[NODES];
static Record   *Ledger[NODES];
static int       ledg_cnt[NODES];
static double    NodePos[NODES][DIM], NodeVel[NODES][DIM], NodeForce[NODES][DIM];
static double    alpha[NODES], beta[NODES], DT_arr[NODES];
static int       adj[NODES][NODES];  // dynamic adjacency

// Messaging queue
#define QSIZE 1024
static QMessage qbuf[QSIZE]; static int qhead, qtail;

// Prototypes
static double hardware_entropy(void);
static double fDyson(double d);
static double wCond(int cnt);
static void   broadcast(int src, Proposal p);
static void   deliver(void);
static void   consensus_step(int node, int step);
static void   compute_forces(void);
static void   integrate_step(void);
static void   log_metrics(int step);
static void   log_tape(int step,int node,int prev_state,int new_state,double ent,double ts);
static void   log_sensors(int step);
static void   rewire_topology(int step);
static void   apply_boundaries(int node);

// Hardware entropy: [0,1)
static double hardware_entropy(void) {
    uint64_t rnd;
    if(getrandom(&rnd,sizeof rnd,0)!=(ssize_t)sizeof rnd){ perror("getrandom"); exit(1);}    
    return (double)rnd/ (double)UINT64_MAX;
}

// Dyson attenuation
static double fDyson(double d) {
    return exp(-Lambda * h_const/(kB_const*T_const*(d+1e-12)));
}

// Occupancy factor
static double wCond(int cnt) {
    double x = cnt % ROUNDS + 1;
    return 1.0/(exp(h_const*x/(kB_const*T_const)) - 1.0);
}

// Enqueue proposal with dynamic adjacency
static void broadcast(int src, Proposal p) {
    for(int dst=0; dst<NODES; dst++) {
        if(dst==src || !adj[src][dst]) continue;
        qbuf[qtail++%QSIZE] = (QMessage){src,dst,p};
    }
}

// Deliver pending
static void deliver(void) {
    while(qhead<qtail) {
        QMessage m = qbuf[qhead++%QSIZE];
        Pending[m.dst][pend_cnt[m.dst]++] = m.p;
    }
}

// Consensus & tape logging
static void consensus_step(int node, int step) {
    deliver();
    if(pend_cnt[node]>0) {
        Proposal pr = Pending[node][0];
        int prev = ledg_cnt[node]>0 ? Ledger[node][ledg_cnt[node]-1].mode : -1;
        Record rec = {pr.mode, ledg_cnt[node]+1, pr.ts, (double)prev, 0.0};
        uint64_t hval = ((uint64_t)(pr.ts*1e6)) ^ (uint64_t)prev;
        rec.hash=(double)hval;
        Ledger[node][ledg_cnt[node]++]=rec;
        log_tape(step,node,prev,rec.mode,pr.ts,step*DT_arr[node]);
    }
    pend_cnt[node]=0;
}

// Compute forces projected to 2D
static void compute_forces(void) {
    for(int i=0;i<NODES;i++) for(int k=0;k<DIM;k++) NodeForce[i][k]=0;
    for(int i=0;i<NODES;i++){
        for(int j=i+1;j<NODES;j++){
            double dx[DIM],dist2=0;
            for(int k=0;k<DIM;k++){ dx[k]=NodePos[j][k]-NodePos[i][k]; dist2+=dx[k]*dx[k]; }
            double d=sqrt(dist2)+1e-12;
            double rep=1.0/(d*d), att=d*d;
            double net=(att-rep)*fDyson(d)*wCond(ledg_cnt[i]+ledg_cnt[j]);
            for(int k=0;k<2;k++){
                double f=net*(dx[k]/d);
                NodeForce[i][k]+=f; NodeForce[j][k]-=f;
            }
        }
    }
}

// Integrate with per-node DT and boundary bounce
static void integrate_step(void) {
    for(int i=0;i<NODES;i++){
        for(int k=0;k<2;k++){
            NodeVel[i][k]+=DT_arr[i]*NodeForce[i][k];
            NodePos[i][k]+=DT_arr[i]*NodeVel[i][k];
        }
        apply_boundaries(i);
    }
}

// Boundary conditions: bounce in [-1,1]
static void apply_boundaries(int i) {
    for(int k=0;k<2;k++){
        if(NodePos[i][k]>BOUNDARY){ NodePos[i][k]=BOUNDARY; NodeVel[i][k]*=-1; }
        if(NodePos[i][k]<-BOUNDARY){ NodePos[i][k]=-BOUNDARY; NodeVel[i][k]*=-1; }
    }
}

// Time-varying topology rewiring
static void rewire_topology(int step) {
    if(step%REWIRE_INTERVAL!=0) return;
    for(int i=0;i<NODES;i++){
        for(int j=0;j<NODES;j++){
            if(i==j) { adj[i][j]=0; continue; }
            adj[i][j]=(rand()%(NODES/2))<1; // ~50% sparse
        }
    }
}

// Log overall metrics
static void log_metrics(int step) {
    int commits=0; for(int i=0;i<NODES;i++) commits+=ledg_cnt[i];
    double sh=0, li=0;
    for(int i=0;i<NODES;i++){ double p=1.0/NODES; sh-=p*log(p);}    
    for(int i=0;i<NODES;i++) li+=NodePos[i][0]*NodeVel[i][1]-NodePos[i][1]*NodeVel[i][0];
    double et=T_const+hardware_entropy()*kB_const;
    FILE*f=fopen("chunk_data.csv","a");
    if(step==1) fprintf(f,"step,commits,shannon,liouville,eff_temp\n");
    fprintf(f,"%d,%d,%.6e,%.6e,%.6f\n",step,commits,sh,li,et);
    fclose(f);
}

// Log tape events
static void log_tape(int step,int node,int prev_state,int new_state,double ent,double ts) {
    FILE*f=fopen("tape_data.csv","a");
    if(step==1&&node==0) fprintf(f,"step,node,prev,new,entropy,time\n");
    fprintf(f,"%d,%d,%d,%d,%.6e,%.6e\n",step,node,prev_state,new_state,ent,ts);
    fclose(f);
}

// Log sensor heat
static void log_sensors(int step) {
    FILE*f=fopen("sensor_data.csv","a");
    if(step==1) fprintf(f,"step,i,j,heat\n");
    for(int i=0;i<NODES;i++){
        for(int j=i+1;j<NODES;j++){
            double heat=fabs(NodeForce[i][0]-NodeForce[j][0])+
                        fabs(NodeForce[i][1]-NodeForce[j][1]);
            fprintf(f,"%d,%d,%d,%.6e\n",step,i,j,heat);
        }
    }
    fclose(f);
}

// Main
int main(void) {
    remove("chunk_data.csv"); remove("tape_data.csv"); remove("sensor_data.csv");
    srand((unsigned)time(NULL));
    // init per-node params and adjacency
    for(int i=0;i<NODES;i++){
        Pending[i]=calloc(MODES,sizeof*Pending[i]); Ledger[i]=calloc(STEPS,sizeof*Ledger[i]);
        pend_cnt[i]=ledg_cnt[i]=0;
        alpha[i]=0.5+0.5*(double)i/NODES;
        beta[i]=0.1+0.05*(double)i/NODES;
        DT_arr[i]=BASE_DT*(1+0.2*((double)rand()/RAND_MAX-0.5));
        for(int j=0;j<NODES;j++) adj[i][j]=(i!=j);
        for(int k=0;k<DIM;k++){ NodePos[i][k]=2.0*((double)rand()/RAND_MAX)-1.0; NodeVel[i][k]=0; }
    }
    for(int step=1;step<=STEPS;step++){
        rewire_topology(step);
        for(int i=0;i<NODES;i++){
            Proposal p={rand()%MODES,hardware_entropy()};
            broadcast(i,p);
        }
        for(int i=0;i<NODES;i++) consensus_step(i,step);
        compute_forces(); integrate_step();
        log_metrics(step); log_sensors(step);
    }
    return 0;
}