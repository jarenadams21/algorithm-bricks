#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>

#ifdef __APPLE__
#include <stdlib.h>   // arc4random_buf
extern void arc4random_buf(void *buf, size_t nbytes);
#else
#include <sys/random.h> // getrandom
#endif

// Physical constants and parameters
#define kB            1.380649e-23    // Boltzmann constant
#define h_P           6.62607015e-34  // Planck constant
#define T_NET0        300.0           // Base network temperature (K)
#define DT            10e-19            // Time step (s)
#define STEPS         1000           // Number of simulation steps
#define NODES         3          // Number of nodes
#define MASS          9.109382e-13    // Mass of each node (kg)
#define THERMAL_SCALE 0.1            // Thermal noise scale
#define J_CONST       100.0            // Entanglement strength
#define ALPHA         1.0             // Attack decay constant
#define GRID_MIN     -10.5            // Grid sampling bounds
#define GRID_MAX      10.5
#define GRID_STEP     0.667             // Coarser grid
#define MSG_COUNT     10000           // Number of message attempts per benchmark

typedef struct { double x, y, vx, vy; } Node;
typedef struct { double x, y; } Point;

//======================================================================
// Hardware entropy → uniform [0,1)
static void get_entropy(void *buf, size_t len) {
#ifdef __APPLE__
    arc4random_buf(buf, len);
#else
    if (getrandom(buf, len, 0) != (ssize_t)len) {
        perror("getrandom");
        exit(EXIT_FAILURE);
    }
#endif
}
static double rand_uniform() {
    uint64_t r;
    get_entropy(&r, sizeof r);
    return (double)r / (double)UINT64_MAX;
}

// Quantum-inspired force
static double quantum_force(double d) {
    if (d < 1e-6) return 0.0;
    double rep   = d*d - 1.0/(d*d);
    double decay = exp(-h_P/(kB*T_NET0*d));
    return rep * decay;
}

// Route probability (Bose-Einstein)
static double route_probability(double E) {
    double e = exp(E/(kB*T_NET0));
    return 1.0/(e - 1.0);
}

// Entanglement reliability
static double entanglement_reliability() {
    return 1.0 / (1.0 + exp(-2.0 * J_CONST / (kB * T_NET0)));
}

// Convex hull (Monotone Chain)
static int cmp_point(const void *a, const void *b) {
    const Point *p = a, *q = b;
    if (p->x < q->x) return -1;
    if (p->x > q->x) return  1;
    if (p->y < q->y) return -1;
    if (p->y > q->y) return  1;
    return 0;
}
static double cross(const Point O, const Point A, const Point B) {
    return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
}
static void convex_hull(Point *pts, int n, Point *H, int *h) {
    qsort(pts, n, sizeof *pts, cmp_point);
    int m = 0;
    // lower hull
    for (int i = 0; i < n; i++) {
        while (m >= 2 && cross(H[m-2], H[m-1], pts[i]) <= 0) m--;
        H[m++] = pts[i];
    }
    // upper hull
    for (int i = n-2, t = m+1; i >= 0; i--) {
        while (m >= t && cross(H[m-2], H[m-1], pts[i]) <= 0) m--;
        H[m++] = pts[i];
    }
    *h = m - 1;
}

// Point-in-polygon (ray casting)
static int point_in_poly(const Point *poly, int n, Point p) {
    int c = 0;
    for (int i = 0, j = n-1; i < n; j = i++) {
        if (((poly[i].y > p.y) != (poly[j].y > p.y)) &&
            (p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y)
                 / (poly[j].y - poly[i].y) + poly[i].x))
        {
            c = !c;
        }
    }
    return c;
}

// Distance from P to segment AB
static double point_segment_dist(Point P, Point A, Point B) {
    double dx = B.x - A.x, dy = B.y - A.y;
    double t  = ((P.x - A.x)*dx + (P.y - A.y)*dy) / (dx*dx + dy*dy);
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    double projx = A.x + t*dx, projy = A.y + t*dy;
    double ddx   = P.x - projx, ddy   = P.y - projy;
    return sqrt(ddx*ddx + ddy*ddy);
}

//======================================================================
int main(void) {
    // --- Initialize nodes
    Node nodes[NODES] = {
        {  1.0,   0.0, 0.0, 0.0 },
        { -0.5,  sqrt(3)/2, 0.0, 0.0 },
        { -0.5, -sqrt(3)/2, 0.0, 0.0 }
    };
    // initial thermal kick
    for (int i = 0; i < NODES; i++) {
        double u1 = rand_uniform(), u2 = rand_uniform();
        double mag = THERMAL_SCALE * sqrt(-2.0 * log(u1));
        nodes[i].vx = mag * cos(2*M_PI*u2);
        nodes[i].vy = mag * sin(2*M_PI*u2);
    }

    // --- Open summary & benchmark CSVs
    FILE *f_sum   = fopen("summary.csv", "w");
    FILE *f_bench = fopen("benchmark.csv", "w");
    if (!f_sum || !f_bench) { perror("fopen"); exit(EXIT_FAILURE); }
    fprintf(f_sum,   "step,P_mean,S_mean_inside,S_mean_outside\n");
    fprintf(f_bench, "step,msg_count,success_in,percent_in,success_out,percent_out,success_cross,percent_cross\n");

    FILE *f_atk = NULL;
    double fx[NODES], fy[NODES];

    // --- Main simulation
    for (int step = 0; step < STEPS; ++step) {
        // zero forces
        for (int i = 0; i < NODES; i++) {
            fx[i] = fy[i] = 0.0;
        }

        // pairwise quantum forces
        for (int i = 0; i < NODES; i++) {
            for (int j = i+1; j < NODES; j++) {
                double dx = nodes[j].x - nodes[i].x;
                double dy = nodes[j].y - nodes[i].y;
                double d  = sqrt(dx*dx + dy*dy);
                double F  = quantum_force(d);
                if (d > 1e-6) {
                    fx[i] += (dx/d)*F;
                    fy[i] += (dy/d)*F;
                    fx[j] -= (dx/d)*F;
                    fy[j] -= (dy/d)*F;
                }
            }
        }

        // integrate motion
        for (int i = 0; i < NODES; i++) {
            nodes[i].vx += fx[i]/MASS * DT;
            nodes[i].vy += fy[i]/MASS * DT;
            nodes[i].x  += nodes[i].vx * DT;
            nodes[i].y  += nodes[i].vy * DT;
        }

        // compute route probabilities & entanglement
        double d01 = hypot(nodes[0].x - nodes[1].x, nodes[0].y - nodes[1].y);
        double d12 = hypot(nodes[1].x - nodes[2].x, nodes[1].y - nodes[2].y);
        double d20 = hypot(nodes[2].x - nodes[0].x, nodes[2].y - nodes[0].y);

        double P01    = route_probability(d01);
        double P12    = route_probability(d12);
        double P20    = route_probability(d20);
        double R      = entanglement_reliability();
        double P_mean = (P01 + P12 + P20) / 3.0;

        // build convex hull of the three nodes
        Point pts[NODES], hull[2*NODES];
        for (int i = 0; i < NODES; i++) {
            pts[i].x = nodes[i].x;
            pts[i].y = nodes[i].y;
        }
        int hn;
        convex_hull(pts, NODES, hull, &hn);

        // aggregate inside/outside attack means
        double sum_in=0, sum_out=0;
        int cnt_in=0, cnt_out=0;
        for (double px = GRID_MIN; px <= GRID_MAX; px += GRID_STEP) {
            for (double py = GRID_MIN; py <= GRID_MAX; py += GRID_STEP) {
                Point P = {px, py};
                // distance to mesh boundary
                double min_d = INFINITY;
                for (int i = 0; i < hn; i++) {
                    int j = (i+1) % hn;
                    double dseg = point_segment_dist(P, hull[i], hull[j]);
                    if (dseg < min_d) min_d = dseg;
                }
                double S_pt = exp(-ALPHA * min_d) * (1.0 - R);
                int inside = point_in_poly(hull, hn, P);
                if (inside) { sum_in += S_pt; cnt_in++; }
                else        { sum_out += S_pt; cnt_out++; }
            }
        }
        double S_mean_in  = cnt_in  ? sum_in / cnt_in   : 0.0;
        double S_mean_out = cnt_out ? sum_out / cnt_out : 0.0;

        // write summary
        fprintf(f_sum, "%d,%.6e,%.6e,%.6e\n",
                step, P_mean, S_mean_in, S_mean_out);

        // benchmark superfluid channels
        int succ_in=0, succ_out=0;
        for (int m = 0; m < MSG_COUNT; m++) {
            if (rand_uniform() < R)                succ_in++;
            if (rand_uniform() < (1.0 - S_mean_out)) succ_out++;
        }
        double pct_in  = 100.0 * succ_in / MSG_COUNT;
        double pct_out = 100.0 * succ_out / MSG_COUNT;
        // cross-zone always fails → 0
        fprintf(f_bench, "%d,%d,%d,%.2f,%d,%.2f,0,0.00\n",
                step, MSG_COUNT,
                succ_in, pct_in,
                succ_out, pct_out);

        // on the final step only, dump attack.csv
        if (step == STEPS - 1) {
            f_atk = fopen("attack.csv", "w");
            fprintf(f_atk, "step,px,py,link,R,S,inside\n");
            for (double px = GRID_MIN; px <= GRID_MAX; px += GRID_STEP) {
                for (double py = GRID_MIN; py <= GRID_MAX; py += GRID_STEP) {
                    Point P = {px, py};
                    double min_d = INFINITY;
                    for (int i = 0; i < hn; i++) {
                        int j = (i+1) % hn;
                        double dseg = point_segment_dist(P, hull[i], hull[j]);
                        if (dseg < min_d) min_d = dseg;
                    }
                    double S_pt = exp(-ALPHA * min_d) * (1.0 - R);
                    int inside = point_in_poly(hull, hn, P);
                    fprintf(f_atk, "%d,%.3f,%.3f,01,%.6e,%.6e,%d\n",
                            step, px, py, R, S_pt, inside);
                    fprintf(f_atk, "%d,%.3f,%.3f,12,%.6e,%.6e,%d\n",
                            step, px, py, R, S_pt, inside);
                    fprintf(f_atk, "%d,%.3f,%.3f,20,%.6e,%.6e,%d\n",
                            step, px, py, R, S_pt, inside);
                }
            }
            fclose(f_atk);
        }
    }

    fclose(f_sum);
    fclose(f_bench);
    return 0;
}
