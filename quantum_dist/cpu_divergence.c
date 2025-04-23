// cpu_divergence.c
// Monitors per-core CPU utilization, regularizes, and computes KL divergence against uniform.
//
// Compile: clang -Wall -O2 cpu_divergence.c -o cpu_divergence
// Run: ./cpu_divergence

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <mach/processor_info.h>
#include <mach/mach_host.h>

#define SLEEP_USEC 1000000   // 1 second
#define EPS        0.05      // regularization weight ε

// Fetch per-core total and idle ticks
static boolean_t get_cpu_ticks(long **out_ticks, natural_t *out_count) {
    host_name_port_t host = mach_host_self();
    processor_info_array_t cpu_info;
    mach_msg_type_number_t count;
    kern_return_t kr = host_processor_info(
        host, PROCESSOR_CPU_LOAD_INFO, out_count,
        &cpu_info, &count
    );
    if (kr != KERN_SUCCESS) return FALSE;

    *out_ticks = malloc(sizeof(long) * count);
    if (!*out_ticks) { vm_deallocate(mach_task_self(), (vm_address_t)cpu_info, count * sizeof(integer_t)); return FALSE; }

    // cpu_info is an array of cpu_count blocks of CPU_STATE_*
    for (natural_t i = 0; i < *out_count; i++) {
        integer_t *base = cpu_info + (i * CPU_STATE_MAX);
        // total = user + sys + nice + idle
        (*out_ticks)[i] = base[CPU_STATE_USER] + base[CPU_STATE_SYSTEM] + base[CPU_STATE_NICE] + base[CPU_STATE_IDLE];
    }
    // stash idle separately by encoding negative idle, so later we subtract
    for (natural_t i = 0; i < *out_count; i++) {
        integer_t *base = cpu_info + (i * CPU_STATE_MAX);
        (*out_ticks)[i] = (*out_ticks)[i] - base[CPU_STATE_IDLE];  // now it's (active ticks)
    }
    vm_deallocate(mach_task_self(), (vm_address_t)cpu_info, count * sizeof(integer_t));
    return TRUE;
}

int main() {
    // 1) Get number of cores d
    int d = 0;
    size_t len = sizeof(d);
    sysctlbyname("hw.ncpu", &d, &len, NULL, 0);
    if (d <= 0) d = 1;

    // 2) Allocate buffers for tick snapshots
    long *prev_ticks = NULL, *curr_ticks = NULL;
    natural_t core_count = 0;

    // Initial snapshot
    if (!get_cpu_ticks(&prev_ticks, &core_count) || core_count != (natural_t)d) {
        fprintf(stderr, "Failed to get initial CPU ticks\n");
        return 1;
    }

    // Precompute expected (uniform) distribution:
    double *p_exp = malloc(sizeof(double)*d);
    double base = 1.0 / d;
    for (int i = 0; i < d; i++) {
        // regularize the uniform too
        p_exp[i] = (1.0 - EPS) * base + EPS * base;
    }

    printf("Sampling every %.1f s  (ε=%.3f, cores=%d)\n\n", SLEEP_USEC/1e6, EPS, d);
    printf("Time\tKL-Divergence\n");
    fflush(stdout);

    // Loop: sample every second, compute divergence
    for (int iter = 1; ; iter++) {
        usleep(SLEEP_USEC);

        // get next snapshot
        if (!get_cpu_ticks(&curr_ticks, &core_count)) {
            fprintf(stderr, "Failed to get CPU ticks\n");
            break;
        }

        // Compute active ticks delta per core
        double sum_active = 0.0;
        double *p_meas = malloc(sizeof(double)*d);
        for (int i = 0; i < d; i++) {
            long delta = curr_ticks[i] - prev_ticks[i];
            if (delta < 0) delta = 0;
            p_meas[i] = (double)delta;
            sum_active += p_meas[i];
        }
        // swap buffers for next iteration
        free(prev_ticks);
        prev_ticks = curr_ticks;
        curr_ticks = NULL;

        // turn into probability distribution (avoid div by zero)
        if (sum_active <= 0) sum_active = 1;
        for (int i = 0; i < d; i++) {
            p_meas[i] /= sum_active;
            // regularize: p' = (1-ε)p + ε/d
            p_meas[i] = (1.0 - EPS) * p_meas[i] + EPS * (1.0/d);
        }

        // compute KL(p_meas' || p_exp)
        double kl = 0.0;
        for (int i = 0; i < d; i++) {
            double a = p_meas[i], b = p_exp[i];
            kl += a * log(a / b);
        }

        printf("%3d\t%.6f\n", iter, kl);
        fflush(stdout);

        free(p_meas);
    }

    free(prev_ticks);
    return 0;
}