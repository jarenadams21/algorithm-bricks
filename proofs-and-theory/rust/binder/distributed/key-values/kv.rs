use rand::Rng;
use rand_distr::{Distribution, Laplace};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

/// Key and value stored in the K-V store
/// In Femur, the server side is "public," but the client's queries must be hidden.
pub type Key = u64;
pub type Value = u64;

/// We store the entire dataset sorted by key. In an actual system, you'd want
/// to keep data in a real DB/LSM tree/disk-based structure, or even an in-memory
/// system like Redis.
#[derive(Clone, Debug)]
pub struct PublicDataset {
    /// The keys are sorted. We can store them along with values in a Vec.
    pub data: Vec<(Key, Value)>,
}

/// Simplified server. It has a `PublicDataset` plus any indexing info needed for
/// the "variable-range PIR" scheme. We also store a single "public-key" if needed.
#[derive(Clone)]
pub struct FemurServer {
    pub dataset: Arc<Mutex<PublicDataset>>,
    // Potentially more fields if you want to store FHE keys, etc.
}

/// A PGM-index node for demonstration. Real code might store slope/intercept, etc.
#[derive(Clone, Debug)]
pub struct PgmLayerModel {
    pub slope: f64,
    pub intercept: f64,
    pub min_key: Key,
    pub max_key: Key,
}

/// PGM-index structure. This is what we store client-side (in real Femur).
#[derive(Clone, Debug)]
pub struct PgmIndex {
    /// We show multiple layers with naive representation.
    pub layers: Vec<Vec<PgmLayerModel>>,
    /// Maximum allowable prediction error for leaf nodes
    pub error_leaf: usize,
    /// The number of keys in the dataset
    pub total_keys: usize,
}

/// The user picks a "distance-based" privacy parameter `dp` and a maximum index distance `D`.
/// Then the noise mechanism ensures `dp-dist-indistinguishability` for any pair of queries within distance D.
#[derive(Copy, Clone, Debug)]
pub struct PrivacyParams {
    pub dp: f64, // e.g. 2^-6 if you want "pretty tight" privacy
    pub dist: usize,
}

/// The user can choose how they want to retrieve data from the server:
/// 1) Plaintext download, or
/// 2) Variable-range PIR
#[derive(Copy, Clone, Debug)]
pub enum RetrieveScheme {
    PlaintextDownload,
    VarRangePIR,
}

/// A query from the client to the server, specifying an obfuscated range
/// plus the chosen retrieval scheme. If using VarRangePIR, we'd also
/// have an encrypted "selector vector" or something like that.
#[derive(Clone, Debug)]
pub struct FemurQuery {
    pub left: usize,
    pub right: usize,
    pub scheme: RetrieveScheme,
    // In a real system, we store the FHE-encrypted "one-hot vector" here:
    // pub encrypted_selector: Vec<u8>,
    // ... or other FHE parameters
}

/// The server returns a "FemurResult". For plaintext downloads, we just store the data vector directly.
/// For VarRangePIR, we might store an encrypted result. Here we show a mock encrypted blob.
#[derive(Clone, Debug)]
pub enum FemurResult {
    PlainResult(Vec<(Key, Value)>),
    MockEncryptedResult(Vec<u8>),
}

// -----------------------------------------------------------------------------
//  Simple PGM-Index Construction
// -----------------------------------------------------------------------------

impl PgmIndex {
    /// Build a trivial "linear model" style PGM with a single layer.
    /// For a real PGM-index, see https://github.com/gvinciguerra/PGM-index
    pub fn build_index(sorted_data: &[(Key, Value)], error_leaf: usize) -> Self {
        // We'll just store the total # of keys and not do a real layered approach for the example.
        // Real PGM algorithm is more sophisticated: scanning from the bottom, building segments, etc.
        let total_keys = sorted_data.len();

        // Create a single dummy layer that maps min..max with slope=1, intercept=0 for demonstration
        // This is obviously not a real PGM. We'll rely on a small error range at search time.
        let min_key = sorted_data.first().map(|(k, _)| *k).unwrap_or(0);
        let max_key = sorted_data.last().map(|(k, _)| *k).unwrap_or(0);

        let model = PgmLayerModel {
            slope: 1.0,
            intercept: 0.0,
            min_key,
            max_key,
        };

        PgmIndex {
            layers: vec![vec![model]],
            error_leaf,
            total_keys,
        }
    }

    /// Predict the *position* (rank) of a key in the sorted array. We return an approximate rank.
    /// For demonstration, we do a naive rank transform: 
    ///   rank = (key - min_key) / (max_key - min_key + 1) * total_keys
    /// Then we clamp to [0, total_keys-1].
    pub fn predict_rank(&self, key: Key) -> usize {
        let layer0 = &self.layers[0][0]; // single model
        if layer0.max_key <= layer0.min_key {
            return 0; // degenerate
        }
        let span = (layer0.max_key - layer0.min_key + 1) as f64;
        let rel = (key.saturating_sub(layer0.min_key)) as f64 / span;
        let pred_f = rel * (self.total_keys as f64);
        let pred = if pred_f < 0.0 {
            0
        } else if pred_f > self.total_keys as f64 {
            self.total_keys - 1
        } else {
            pred_f.round() as usize
        };
        pred
    }
}

// -----------------------------------------------------------------------------
//  Noise Generation (Discrete Laplace) for distance-based indistinguishability
// -----------------------------------------------------------------------------

/// Sample from a discrete Laplace distribution with scale = scale.
/// This is a naive approach.  For distance-based privacy, see the paper
/// for details. We do "modular" so negative noise wraps in a certain range.
fn sample_discrete_laplace<R: Rng>(rng: &mut R, scale: f64) -> i64 {
    // We use the continuous Laplace from rand_distr, then round.
    // Real discrete-laplace sampling can be done by custom code, but this is a simpler approximation.
    let laplace = Laplace::new(0.0, scale).unwrap();
    let sample = laplace.sample(rng);
    sample.round() as i64
}

/// Generate the obfuscated range [left, right], containing [pred-error_leaf, pred+error_leaf].
/// We show a naive approach to use discrete-laplace noise and clamp to [0..total_len-1].
pub fn generate_obfuscated_range(
    predicted: usize,
    error_leaf: usize,
    privacy: PrivacyParams,
    total_len: usize,
) -> (usize, usize) {
    // The core region is [pred-error_leaf, pred+error_leaf].
    let base_left = predicted.saturating_sub(error_leaf);
    let base_right = (predicted + error_leaf).min(total_len.saturating_sub(1));

    // We sample noise for left boundary and right boundary with scale = 2 * dist / dp
    // as discussed in the paper. The smaller dp => bigger noise. Larger dist => bigger noise.
    let scale = 2.0 * (privacy.dist as f64) / privacy.dp;

    let mut rng = rand::thread_rng();
    let noise_left = sample_discrete_laplace(&mut rng, scale);
    let noise_right = sample_discrete_laplace(&mut rng, scale);

    // We add the noise to the base_left and base_right. Then clamp.
    // For "distance-based" we want the final guaranteed coverage. 
    // In practice: left = base_left - noise_left, right = base_right + noise_right.
    // But do saturating cast to avoid negative. Then clamp to total_len.
    let mut new_left = if noise_left.is_negative() {
        base_left.saturating_sub(noise_left.unsigned_abs() as usize)
    } else {
        base_left.saturating_add(noise_left as usize)
    };
    let mut new_right = if noise_right.is_negative() {
        base_right.saturating_sub(noise_right.unsigned_abs() as usize)
    } else {
        base_right.saturating_add(noise_right as usize)
    };

    // clamp
    if new_left >= total_len {
        new_left = total_len.saturating_sub(1);
    }
    if new_right >= total_len {
        new_right = total_len.saturating_sub(1);
    }
    if new_right < new_left {
        // swap or unify them so there's at least one valid element
        std::mem::swap(&mut new_left, &mut new_right);
    }

    (new_left, new_right)
}

// -----------------------------------------------------------------------------
//  Server-Side Query Processing
// -----------------------------------------------------------------------------

impl FemurServer {
    pub fn new(dataset: Vec<(Key, Value)>) -> Self {
        let ds = PublicDataset { data: dataset };
        FemurServer {
            dataset: Arc::new(Mutex::new(ds)),
        }
    }

    /// Perform the query according to the given scheme. The server does *not*
    /// learn the real key. It only sees the obfuscated [L,R].
    pub fn handle_query(&self, query: FemurQuery) -> FemurResult {
        match query.scheme {
            RetrieveScheme::PlaintextDownload => self.handle_plaintext(query.left, query.right),
            RetrieveScheme::VarRangePIR => self.handle_varrange_pir(query.left, query.right),
        }
    }

    /// For plaintext download, just slice from [left..=right], and return it.
    fn handle_plaintext(&self, left: usize, right: usize) -> FemurResult {
        let locked = self.dataset.lock().unwrap();
        let total_len = locked.data.len();
        if total_len == 0 {
            return FemurResult::PlainResult(vec![]);
        }

        let l = left.min(total_len - 1);
        let r = right.min(total_len - 1);
        let subset = locked.data[l..=r].to_vec();
        FemurResult::PlainResult(subset)
    }

    /// A mocked "variable-range PIR". We pretend we do homomorphic computations
    /// over the range [left..=right]. The result is a single ciphertext. Real code
    /// would do a homomorphic inner product or use a specialized library.
    fn handle_varrange_pir(&self, left: usize, right: usize) -> FemurResult {
        let locked = self.dataset.lock().unwrap();
        let total_len = locked.data.len();
        if total_len == 0 {
            return FemurResult::MockEncryptedResult(vec![]);
        }
        let l = left.min(total_len - 1);
        let r = right.min(total_len - 1);

        // We are *pretending* to produce a ciphertext that encloses all pairs in [l..=r].
        // In real code, you'd do a homomorphic encoding. This might be large, but for
        // demonstration let's just produce a mock, e.g. a JSON of the data (still not secure).
        // Then you'd encrypt it or do real PIR. We'll just do a silly "serialized" version.
        let subset = locked.data[l..=r].to_vec();
        let serialized = serde_json::to_vec(&subset).unwrap();
        // In real code: encrypt "serialized" with FHE, or do partial expansions, etc.

        FemurResult::MockEncryptedResult(serialized)
    }
}

// -----------------------------------------------------------------------------
//  Client-Side Femur Implementation (key→position→obfuscated range→retrieve)
// -----------------------------------------------------------------------------

/// The client object. It stores a local copy of the PGM-index and a pointer to the
/// remote server. In real usage, you'd have separate processes, networking, etc.
pub struct FemurClient {
    pub index: PgmIndex,
    pub server: FemurServer,
    pub bandwith_mbps: f64, // used for cost-based scheme selection
    pub compute_time_pir_per_plaintext_ms: f64, // a toy example
}

/// A mini cost-model: choose PlaintextDownload if the expected data size is smaller
/// than the cost of homomorphic computations, else use VarRangePIR.
impl FemurClient {
    pub fn new(index: PgmIndex, server: FemurServer) -> Self {
        FemurClient {
            index,
            server,
            bandwith_mbps: 50.0,
            compute_time_pir_per_plaintext_ms: 300.0,
        }
    }

    /// Decide which scheme to use, given an obfuscated range length. This is a *simple*
    /// cost model. The paper discusses a more precise formula with (Comm / bandwidth + computation).
    pub fn choose_retrieve_scheme(&self, range_len: usize) -> RetrieveScheme {
        // Estimate plain download cost:
        //   comm_time (seconds) = (range_len * 16 bytes) / (bandwidth bits => MBps)
        //   ignoring overhead
        let bytes = (range_len as f64) * 16.0; // each key-value pair 16 bytes, for example
        // bandwidth in Mbps => MBps = bandwith_mbps / 8
        let comm_time_plain_sec = bytes / ((self.bandwith_mbps / 8.0) * 1_000_000.0); 
        let comm_time_plain_ms = comm_time_plain_sec * 1000.0;

        // For VarRangePIR, we pretend the data is always 1 ciphertext.
        // So we pay ~1 MB worth of data, plus server computation time.
        // Let’s say 1 MB -> 8 Megabits, plus some constant overhead. Then add compute_time.
        let comm_time_pir_ms = 8.0 / self.bandwith_mbps * 1000.0;
        let compute_time_pir_ms = self.compute_time_pir_per_plaintext_ms;

        let total_plain = comm_time_plain_ms; 
        let total_pir = comm_time_pir_ms + compute_time_pir_ms;

        if total_plain < total_pir {
            RetrieveScheme::PlaintextDownload
        } else {
            RetrieveScheme::VarRangePIR
        }
    }

    /// Public “lookup” method that the client calls to retrieve the value for `key`
    /// at security level (dp, dist).
    ///
    /// Return either the single matching record, or None if not found.
    pub fn lookup(&self, key: Key, pp: PrivacyParams) -> Option<Value> {
        // 1) Key→Position
        let pred_rank = self.index.predict_rank(key);
        // 2) Expand rank with +/- self.index.error_leaf
        //    => [pred_rank - e, pred_rank + e]
        // 3) Add noise for distance-based indistinguishability
        let (obf_left, obf_right) = generate_obfuscated_range(
            pred_rank,
            self.index.error_leaf,
            pp,
            self.index.total_keys,
        );
        let range_len = obf_right.saturating_sub(obf_left) + 1;

        // 4) Cost-based scheme selection
        let scheme = self.choose_retrieve_scheme(range_len);

        // 5) Construct query
        let query = FemurQuery {
            left: obf_left,
            right: obf_right,
            scheme,
        };

        // 6) Send query to server
        let resp = self.server.handle_query(query);

        // 7) Locally parse the result
        match resp {
            FemurResult::PlainResult(pairs) => {
                // Just look for the key in pairs
                pairs.iter().find_map(|(k,v)| if *k == key { Some(*v) } else { None })
            }
            FemurResult::MockEncryptedResult(blob) => {
                // Real code: decrypt the ciphertext / do local filtering.
                // Our “blob” is just a JSON serialization, so parse that:
                let subset: Vec<(Key, Value)> = serde_json::from_slice(&blob).unwrap_or_default();
                subset.iter().find_map(|(k,v)| if *k == key { Some(*v) } else { None })
            }
        }
    }
}

// -----------------------------------------------------------------------------
//  Example main that constructs a dataset, builds a PGM-index, and does queries
// -----------------------------------------------------------------------------

fn main() {
    // Suppose we have a dataset of 20 key-value pairs, sorted by key
    let mut dataset: Vec<(Key, Value)> = vec![];
    for i in 0..20u64 {
        dataset.push((i * 10, i + 1000));
    }
    // Sort by key, though we already inserted in ascending order
    dataset.sort_by_key(|&(k, _)| k);

    // Build a server
    let server = FemurServer::new(dataset.clone());

    // Build a PGM-index on the client side with an error_leaf=2
    let pgm_index = PgmIndex::build_index(&dataset, 2);

    // Initialize our client
    let client = FemurClient::new(pgm_index, server);

    // Suppose the user wants to query the key=55. The real key’s rank
    // is around position=5 or 6. We pick a "distance"=10, dp=1/64.0 (2^-6).
    let pp = PrivacyParams { dp: 1.0/64.0, dist: 10 };

    // Do the lookup
    let val = client.lookup(55, pp);
    println!("Queried key=55 => Value = {:?}", val);

    // Let’s do another lookup with bigger distance => more noise => bigger obfuscated range
    let pp2 = PrivacyParams { dp: 1.0/64.0, dist: 100 };
    let val2 = client.lookup(55, pp2);
    println!("Queried key=55 with bigger privacy => Value = {:?}", val2);

    // The server does not know we asked for key=55, it only sees an obfuscated range
    // that might be quite large, depending on the noise.
}
