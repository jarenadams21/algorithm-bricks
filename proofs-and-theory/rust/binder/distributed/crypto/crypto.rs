use rand::Rng;
use serde::Deserialize;
use std::fs;

//
// 1) JSON Parsing
//

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
enum AccessNode {
    Party {
        id: String,
    },
    Threshold {
        k: usize,
        children: Vec<AccessNode>,
    },
}

const P: i64 = 97;

fn modp(x: i64) -> i64 {
    let r = x % P;
    if r < 0 { r + P } else { r }
}
fn addp(a: i64, b: i64) -> i64 { modp(a + b) }
fn mulp(a: i64, b: i64) -> i64 { modp(a * b) }

// Extended Euclid for inverse mod p
fn inv_modp(a: i64) -> i64 {
    let (g, x, _) = extended_gcd(a, P);
    if g.abs() != 1 {
        panic!("No inverse for {} mod {}", a, P);
    }
    modp(x)
}
fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        return (a, 1, 0);
    }
    let (g2, x2, y2) = extended_gcd(b, a % b);
    (g2, y2, x2 - (a / b) * y2)
}

fn parse_access_structure(filename: &str) -> AccessNode {
    let data = fs::read_to_string(filename)
        .expect("Unable to read JSON file");
    serde_json::from_str::<AccessNode>(&data)
        .expect("JSON parse error")
}

//
// 2) Monotone Span Program (MSP) via insertion-based
//

#[derive(Debug)]
struct Msp {
    matrix: Vec<Vec<i64>>, // (m x d)
    row_to_party: Vec<String>,
}

/// A standard n×k Vandermonde for x_i = i+1 => row i => [1, x, x^2, ...]
fn build_vandermonde(n: usize, k: usize) -> Vec<Vec<i64>> {
    let mut mat = vec![vec![0; k]; n];
    for (i, row) in mat.iter_mut().enumerate() {
        let x = (i as i64) + 1; // distinct x in {1..n}
        let mut cur = 1;
        for col in 0..k {
            row[col] = modp(cur);
            cur = mulp(cur, x);
        }
    }
    mat
}

/// Build the MSP from a node
fn build_msp(node: &AccessNode) -> Msp {
    match node {
        AccessNode::Party { id } => {
            // Single row => [1]
            Msp {
                matrix: vec![vec![1]],
                row_to_party: vec![id.clone()],
            }
        }
        AccessNode::Threshold { k, children } => {
            insertion_threshold(*k, children)
        }
    }
}

/// insertion_threshold(k, children):
/// 1) Build an n×k Vandermonde (n = #children)
/// 2) Start with that as the parent's matrix => n rows
/// 3) For each child i:
///    - remove row i
///    - expand columns if needed
///    - multiply child's last k columns by row i of Vandermonde
///    - insert child's rows
fn insertion_threshold(k: usize, children: &[AccessNode]) -> Msp {
    let n = children.len();
    // If no children => trivial
    if n == 0 {
        return Msp {
            matrix: vec![],
            row_to_party: vec![],
        };
    }
    // Build child MSPs
    let child_msps: Vec<Msp> = children.iter().map(build_msp).collect();

    // Build parent Vandermonde: n×k
    let vm = build_vandermonde(n, k);

    // Start with parent_msp = Vandermonde
    // row labels => "placeholder"
    let mut parent_matrix = vm.clone();
    let mut parent_labels = vec!["placeholder".to_string(); n];

    let mut parent_msp = Msp {
        matrix: parent_matrix,
        row_to_party: parent_labels,
    };

    // We'll do insertion in ascending order of i
    let mut offset = 0;
    for i in 0..n {
        // row in parent's MSP is i + offset
        let row_index = i + offset;
        let row_vand = &vm[i]; // e.g. [1, 2, 4, ...] for row i

        // child i's MSP
        let c_msp = &child_msps[i];

        // remove that row => effectively "free up" the slot
        parent_msp.matrix.remove(row_index);
        parent_msp.row_to_party.remove(row_index);

        // expand parent's columns if needed
        let parent_cols = if parent_msp.matrix.is_empty() {
            0
        } else {
            parent_msp.matrix[0].len()
        };
        let child_cols = if c_msp.matrix.is_empty() {
            0
        } else {
            c_msp.matrix[0].len()
        };
        let needed = parent_cols.max(child_cols).max(k);

        // pad parent horizontally
        if parent_cols < needed {
            for rowp in parent_msp.matrix.iter_mut() {
                rowp.resize(needed, 0);
            }
        }

        // copy child's matrix => pad horizontally
        let mut c_mat_expanded = c_msp.matrix.clone();
        for row_c in c_mat_expanded.iter_mut() {
            row_c.resize(needed, 0);
        }

        // multiply child's last k columns by row i's Vandermonde
        let start_col = needed - k; // index of the last k columns
        for row_c in c_mat_expanded.iter_mut() {
            for col_k in 0..k {
                let f = row_vand[col_k];
                row_c[start_col + col_k] = mulp(row_c[start_col + col_k], f);
            }
        }

        // insert child's rows => row_index
        for (r, c_row) in c_mat_expanded.into_iter().enumerate() {
            parent_msp.matrix.insert(row_index + r, c_row);
        }
        for (r, lbl) in c_msp.row_to_party.iter().enumerate() {
            parent_msp.row_to_party.insert(row_index + r, lbl.clone());
        }

        // net gain => child_msp.rows - 1
        offset += c_msp.matrix.len() - 1;
    }

    parent_msp
}

//
// 3) Secret Sharing & Reconstruction
//

fn share_secret(msp: &Msp, secret: i64) -> (Vec<i64>, Vec<i64>) {
    let m = msp.matrix.len();
    let d = if m == 0 { 1 } else { msp.matrix[0].len() };

    // r = (secret, random, ...)
    let mut rng = rand::rng();
    let mut coeff = vec![0i64; d];
    coeff[0] = modp(secret);
    for i in 1..d {
        coeff[i] = modp(rng.random::<i64>().abs());
    }

    let mut shares = vec![0i64; m];
    for i in 0..m {
        let mut s = 0;
        for j in 0..d {
            s = addp(s, mulp(msp.matrix[i][j], coeff[j]));
        }
        shares[i] = s;
    }
    (shares, coeff)
}

fn reconstruct_secret(msp: &Msp, subset_rows: &[usize], subset_shares: &[i64]) -> i64 {
    if msp.matrix.is_empty() {
        return 0;
    }
    let d = msp.matrix[0].len();
    let a_len = subset_rows.len();

    let mut m_a = vec![vec![0i64; d]; a_len];
    for (i, &ridx) in subset_rows.iter().enumerate() {
        m_a[i] = msp.matrix[ridx].clone();
    }

    // build m_a^T
    let mut m_a_t = vec![vec![0i64; a_len]; d];
    for i in 0..a_len {
        for j in 0..d {
            m_a_t[j][i] = m_a[i][j];
        }
    }

    // solve m_a_t * lambda = e1
    for i in 0..d {
        m_a_t[i].push(if i == 0 { 1 } else { 0 });
    }
    let cols = a_len + 1;

    let mut row = 0;
    // forward-elim
    for col in 0..a_len {
        if row >= d {
            break;
        }
        let mut pivot = row;
        while pivot < d && m_a_t[pivot][col] == 0 {
            pivot += 1;
        }
        if pivot == d {
            continue;
        }
        if pivot != row {
            m_a_t.swap(row, pivot);
        }
        let invp = inv_modp(m_a_t[row][col]);
        for c2 in col..cols {
            m_a_t[row][c2] = mulp(m_a_t[row][c2], invp);
        }
        for r2 in (row+1)..d {
            let f = m_a_t[r2][col];
            if f != 0 {
                for c2 in col..cols {
                    let x = m_a_t[r2][c2] - mulp(m_a_t[row][c2], f);
                    m_a_t[r2][c2] = modp(x);
                }
            }
        }
        row += 1;
    }

    // back-sub
    for r in (0..row).rev() {
        let mut pivot_col = 0;
        while pivot_col < a_len && m_a_t[r][pivot_col] == 0 {
            pivot_col += 1;
        }
        if pivot_col == a_len {
            continue;
        }
        for r2 in 0..r {
            let f = m_a_t[r2][pivot_col];
            if f != 0 {
                for c2 in pivot_col..cols {
                    let x = m_a_t[r2][c2] - mulp(m_a_t[r][c2], f);
                    m_a_t[r2][c2] = modp(x);
                }
            }
        }
    }

    // read solution
    let mut lambda = vec![0i64; a_len];
    let mut pivot_row = 0;
    for col in 0..a_len {
        if pivot_row < d && m_a_t[pivot_row][col] == 1 {
            lambda[col] = m_a_t[pivot_row][a_len];
            pivot_row += 1;
        }
    }

    let mut s = 0;
    for (li, &shr) in lambda.iter().zip(subset_shares) {
        s = addp(s, mulp(*li, shr));
    }
    s
}

//
// 4) MAIN
//

fn main() {
    let as_tree = parse_access_structure("access_structure.json");
    println!("Parsed Access Tree: {:#?}", as_tree);

    let msp = build_msp(&as_tree);

    println!("\n==== MSP ====\nMatrix M (mod {}):", P);
    for (i, row) in msp.matrix.iter().enumerate() {
        println!(" Row {i} => party {}: {:?}", msp.row_to_party[i], row);
    }

    let secret = 42;
    println!("\n==== SHARE SECRET = {} (mod {}) ====", secret, P);
    let (all_shares, coeff) = share_secret(&msp, secret);
    println!("Random polynomial coeff = {:?}", coeff);
    println!("All shares             = {:?}", all_shares);

    // Let's test {P1,P2} vs. {P1,P2,P3}
    // For "2-of-[ P1, (2-of-[P2,P3,P4]) ]":

    // a) {P1,P2} => child #2 only has P2 => 1-of-3 => not satisfied => should FAIL
    let mut rows_12 = Vec::<usize>::new();
    let mut shares_12 = Vec::<i64>::new();
    for (i, party) in msp.row_to_party.iter().enumerate() {
        if party == "P1" || party == "P2" {
            rows_12.push(i);
            shares_12.push(all_shares[i]);
        }
    }
    println!("\nTrying reconstruction from {{P1,P2}} => rows {rows_12:?}");
    let rec_12 = reconstruct_secret(&msp, &rows_12, &shares_12);
    println!("Reconstructed = {rec_12}");
    if rec_12 == modp(secret) {
        println!("** Should have failed, but got success!?");
    } else {
        println!("Correctly failed => {rec_12} != {secret}");
    }

    // b) {P1,P2,P3} => child #2 has P2,P3 => 2-of-3 => satisfied => entire threshold => success
    let mut rows_123 = Vec::<usize>::new();
    let mut shares_123 = Vec::<i64>::new();
    for (i, party) in msp.row_to_party.iter().enumerate() {
        if party == "P1" || party == "P2" || party == "P3" {
            rows_123.push(i);
            shares_123.push(all_shares[i]);
        }
    }
    println!("\nTrying reconstruction from {{P1,P2,P3}} => rows {rows_123:?}");
    let rec_123 = reconstruct_secret(&msp, &rows_123, &shares_123);
    println!("Reconstructed = {rec_123}");
    if rec_123 == modp(secret) {
        println!("Success => got {secret}!");
    } else {
        println!("** Should have succeeded => mismatch => {rec_123}");
    }
}
