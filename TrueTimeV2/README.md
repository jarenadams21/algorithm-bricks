**Objective:**
Designed a quantum memory tape and vortex-based phase-space overlay that encodes state transitions of the mesh protocol, and visualize it as a swirling (vortex) pattern in the Python visualizer.

---

## 1. Quantum Memory Tape

1. **Tape Structure**
   - A 1D array `Tape[0..L-1]`, where each cell records a state transition event.
   - Each entry: `[step, node_id, prev_state, new_state, entropy, time_stamp]`.
   - `L` is total number of consensus events (i.e., sum of all `ledgCnt` across steps).

2. **Encoding Transitions**
   - On each consensus commit at node `i`, append to the tape:
     - `step` (iteration count)
     - `node_id = i`
     - `prev_state = last mode in Ledger[i] - 1` (or 0 at first)
     - `new_state = current mode`
     - `entropy = hardware_entropy()`
     - `time_stamp = global simulation time (step * DT)`

## 2. Phase-Space Vortex Overlay

1. **Phase Coordinates**
   - Map each tape entry to a 2D point `(x, y)` in phase-space:
     - `x = cos(2π * (entropy mod 1)) * R(step)`
     - `y = sin(2π * (entropy mod 1)) * R(step)`
     - `R(step) = α * sqrt(step)` (spiral radius growth), with α scaling factor.

2. **Vortex Pattern**
   - As `step` increases, points trace out a swirling spiral (vortex).
   - Color-code by `node_id` or `new_state`.
   - Optionally, draw arrows/tubes indicating transition direction (prev -> new).

## 3. Python Visualizer Updates

1. **Load Tape**
   - Extend `load_data()` to read `tape_data.csv` (columns: `step,node,prev_state,new_state,entropy,time_stamp`).

2. **Compute Vortex Coordinates**
   - For each row:
     ```python
     theta = 2 * np.pi * (row.entropy % 1)
     r = alpha * np.sqrt(row.step)
     x = r * np.cos(theta)
     y = r * np.sin(theta)
     ```

3. **Plot Spiral**
   - Scatter plot `(x, y)` points with:
     - `c=node_id`, `cmap='viridis'` or state-based colors.
     - `s=5` for point size.
   - Overlay line segments connecting successive points for same node to show flow.

4. **Integrate with Existing Metrics**
   - Add a subplot in the existing figure window for the vortex panel.
   - Layout: use `plt.subplot(2,2,4)` for vortex.

## 4. Data Export

1. Modify C simulator to write `tape_data.csv` as events occur:
   - After each `ledgCnt[node]++`, write a CSV line for that transition.
   - Header: `step,node,prev_state,new_state,entropy,time_stamp`.

2. Ensure `tape_data.csv` resets at simulation start.

---

Once this spec is approved, we can implement the C exporter and update the Python visualizer accordingly.