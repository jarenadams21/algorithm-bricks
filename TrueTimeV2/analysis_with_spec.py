import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -----------------------------------------------------------------------------
# Phase‐Space, Physical 3D, Spectral & Hardware Visualization for Stakeholders
# -----------------------------------------------------------------------------

def load_data():
    df_metrics   = pd.read_csv('metrics.csv')
    df_hw        = pd.read_csv('hardware_stats.csv')
    df_chunk     = pd.read_csv('chunk_data.csv')[['step', 'eff_temp']]
    df_positions = pd.read_csv('position_data.csv')
    df_spectral  = pd.read_csv('spectral.csv')
    # ensure numeric steps and spectral columns
    for df in (df_metrics, df_chunk, df_positions, df_hw, df_spectral):
        df['step'] = pd.to_numeric(df['step'], errors='coerce').astype(int)
        df.dropna(subset=['step'], inplace=True)
        df['step'] = df['step'].astype(int)
    df_spectral['lambda1'] = pd.to_numeric(df_spectral['lambda1'], errors='coerce')
    df_spectral['lambda2'] = pd.to_numeric(df_spectral['lambda2'], errors='coerce')
    return df_metrics, df_chunk, df_positions, df_hw, df_spectral

def compute_defense_metrics(df):
    df = df.copy()
    p = df['edge_frac'].clip(1e-12, 1-1e-12)
    df['edge_entropy'] = -p * np.log(p) - (1-p) * np.log(1-p)
    baseline = 0.5
    df['strength'] = df['safety'] / baseline
    return df

def plot_physical_snapshots(df_chunk, df_positions, key_steps):
    df = df_positions.merge(df_chunk, on='step', how='left')
    fig = plt.figure(figsize=(15,5))
    for idx, step in enumerate(key_steps, start=1):
        ax = fig.add_subplot(1,3,idx, projection='3d')
        sub = df[df['step']==step]
        sc = ax.scatter(sub['x'], sub['y'], sub['z'],
                        c=sub['eff_temp'], cmap='inferno', s=40)
        ax.set_title(f'Step {step}')
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        fig.colorbar(sc, ax=ax, shrink=0.6, label='Eff Temp')
    plt.suptitle('3D Positions at Collapse / Mid / Superfluid', y=1.02)
    plt.tight_layout()
    plt.show()

def plot_phase_space_evolution(df_metrics, key_steps):
    x = df_metrics['safety']
    y = df_metrics['liveness']
    t = df_metrics['step']
    cmap = plt.get_cmap('viridis')

    fig, ax = plt.subplots(figsize=(6,6))
    for i in range(len(df_metrics)-1):
        c = cmap((t.iat[i]-t.min())/(t.max()-t.min()))
        ax.plot(x.iloc[i:i+2], y.iloc[i:i+2], color=c, linewidth=2)

    for step, marker, color in zip(key_steps, ['o','s','*'], ['red','orange','green']):
        row = df_metrics[df_metrics['step']==step]
        if not row.empty:
            ax.scatter(row['safety'], row['liveness'],
                       s=150, marker=marker, edgecolor='k', color=color,
                       label=f'Step {step}')

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=t.min(), vmax=t.max()))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label='Step')
    ax.set_xlabel('Safety'); ax.set_ylabel('Liveness')
    ax.set_title('Phase‐Space: Safety vs Liveness')
    ax.legend(loc='lower right')
    ax.grid(True)
    plt.tight_layout()
    plt.show()

def plot_defense_time_overlap(df_defense, key_steps):
    time = df_defense['step']
    plt.figure(figsize=(10,4))
    plt.plot(time, df_defense['safety'],    label='Safety',    linewidth=1.5)
    plt.plot(time, df_defense['liveness'],  label='Liveness',  linewidth=1.5)
    plt.plot(time, df_defense['edge_frac'], label='Edge Fraction', linewidth=1.5)
    ee_norm = df_defense['edge_entropy'] / df_defense['edge_entropy'].max()
    plt.plot(time, ee_norm,                label='Edge Entropy (norm)', linewidth=1.5)
    plt.plot(time, df_defense['strength'], label='Strength',  linewidth=1.5)

    for step, color in zip(key_steps, ['red','orange','green']):
        plt.axvline(step, color=color, linestyle='--', linewidth=1.5,
                    label=f'Key Step {step}')

    plt.xlabel('Step'); plt.ylabel('Value')
    plt.title('Defense Metrics Over Time with Key Events')
    plt.legend(ncol=3, loc='upper center')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_spectral_gap(df_spectral, key_steps):
    time = df_spectral['step']
    gap  = df_spectral['lambda1'] - df_spectral['lambda2']

    plt.figure(figsize=(10,3))
    plt.plot(time, gap, label='Spectral Gap', linewidth=2)
    for step, color in zip(key_steps, ['red','orange','green']):
        plt.axvline(step, color=color, linestyle='--', linewidth=1.5)
    plt.xlabel('Step'); plt.ylabel('λ₁ - λ₂')
    plt.title('Spectral Gap Over Time')
    plt.legend(['Spectral Gap'] + [f'Step {s}' for s in key_steps])
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_hardware_time_series(df_hw, key_steps):
    plt.figure(figsize=(10,4))
    plt.plot(df_hw['step'], df_hw['cpu_temp_c'],    label='CPU Temp (°C)', linewidth=1)
    plt.plot(df_hw['step'], df_hw['rng_latency_s'], label='RNG Latency (s)', linewidth=1)
    plt.plot(df_hw['step'], df_hw['clock_jitter_s'],label='Clock Jitter (s)', linewidth=1)
    for step, color in zip(key_steps, ['red','orange','green']):
        plt.axvline(step, color=color, linestyle='--', linewidth=1)
    plt.xlabel('Step'); plt.ylabel('Value')
    plt.title('Hardware Metrics with Key Events')
    plt.legend(loc='upper right', ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_hardware_vs_strength(df_defense, df_hw, key_steps):
    df = df_defense.merge(df_hw, on='step')
    plt.figure(figsize=(6,4))
    sc = plt.scatter(df['cpu_temp_c'], df['strength'],
                     c=df['step'], cmap='plasma', s=60)
    for step, marker, color in zip(key_steps, ['o','s','*'], ['red','orange','green']):
        row = df[df['step']==step]
        if not row.empty:
            plt.scatter(row['cpu_temp_c'], row['strength'],
                        s=200, marker=marker, edgecolor='k', color=color,
                        label=f'Step {step}')
    plt.colorbar(sc, label='Step')
    plt.xlabel('CPU Temp (°C)'); plt.ylabel('Strength')
    plt.title('Strength vs CPU Temperature')
    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    df_metrics, df_chunk, df_positions, df_hw, df_spectral = load_data()
    df_defense = compute_defense_metrics(df_metrics)

    max_step        = df_metrics['step'].max()
    collapse_step   = int(0.4 * max_step)
    superfluid_step = int(0.9 * max_step)
    mid_step        = (collapse_step + superfluid_step) // 2
    key_steps       = [collapse_step, mid_step, superfluid_step]

    plot_physical_snapshots(df_chunk, df_positions, key_steps)
    plot_phase_space_evolution(df_defense, key_steps)
    plot_defense_time_overlap(df_defense, key_steps)
    plot_spectral_gap(df_spectral, key_steps)
    plot_hardware_time_series(df_hw, key_steps)
    plot_hardware_vs_strength(df_defense, df_hw, key_steps)
