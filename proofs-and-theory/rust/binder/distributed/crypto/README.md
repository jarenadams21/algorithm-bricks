# Crypto.rs

1. Each party may have a different security profile (e.g., running on different operating systems or in different physical locations). You want certain subsets of parties to be “authorized” to reconstruct a secret or to perform some distributed action, and other subsets to be “unauthorized.”

2. A monotone Boolean formula or a threshold expression captures precisely which subsets of parties are allowed. The Rust code:

``` Parses this JSON file into an internal tree data structure. Builds a Monotone Span Program (MSP) from it. Runs secret-sharing (or other distributed-crypto tasks) according to that MSP. If you do not provide the file, the code panics because it cannot find the structure describing how to distribute shares. ```

3. Standard references on insertion-based MSP constructions:
- C. G. Gao, G. J. Simmons. “On the Complexity of Multiparty Computation with Spies.” (1996)
- D. Beaver, S. Micali, P. Rogaway. “The Round Complexity of Secure Protocols.” STOC (1990)
- M. Karchmer and A. Wigderson. “On Span Programs.” (1989)

## Needs Work
Threshold children requirements are still not being respected by the algorithm:

Trying reconstruction from {P1,P2} => rows [0, 1]
Reconstructed = 42
** Should have failed, but got success!?

Trying reconstruction from {P1,P2,P3} => rows [0, 1, 2]
Reconstructed = 42
Success => got 42!