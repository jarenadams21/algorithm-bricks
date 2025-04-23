/*
 * Atomically safe time to pursue convergence
 */
fn main() {
    // Dimension of the system
    let d = 3usize;
    // Noise weight ε
    let eps: f64 = 0.2;
    // Threshold δ
    let delta: f64 = 0.3;

    // Initial eigenvalues of ρ
    let mut eigen = vec![0.6, 0.3, 0.1];
    let lam0 = *eigen
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    println!("Iteration\tEmpirical λ_min\tTheoretical λ_min");
    let mut k = 0;
    loop {
        // Empirical minimum eigenvalue
        let lam_emp = *eigen
            .iter()
            .fold(&std::f64::INFINITY, |a, b| if b < a { b } else { a });
        // Theoretical minimum via closed-form
        let lam_th = lam0 * (1.0 - eps).powi(k) + (1.0 - (1.0 - eps).powi(k)) / (d as f64);

        println!("{}\t{:.6}\t{:.6}", k, lam_emp, lam_th);

        if lam_emp >= delta {
            println!(
                "\nReached λ_min >= δ = {} at iteration {} (λ_min = {:.6})",
                delta, k, lam_emp
            );
            break;
        }

        // Apply E(ρ) = (1 - ε)ρ + ε (I/d)
        for x in eigen.iter_mut() {
            *x = (1.0 - eps) * (*x) + eps * (1.0 / (d as f64));
        }
        k += 1;
        if k > 100 {
            eprintln!("Exceeded max iterations");
            break;
        }
    }
}

/*
Sample output when running `cargo run`:
Iteration   Empirical λ_min    Theoretical λ_min
0           0.100000           0.100000
1           0.146667           0.146667
2           0.184000           0.184000
3           0.213867           0.213867
4           0.237760           0.237760
5           0.256875           0.256875
6           0.272166           0.272166
7           0.284400           0.284400
8           0.294186           0.294186
9           0.302016           0.302016

Reached λ_min >= δ = 0.3 at iteration 9 (λ_min = 0.302016)
*/
