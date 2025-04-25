import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Quantum Mesh Vortex Metrics Visualizer
# -----------------------------------------------------------------------------
# Reads:
#   - 'chunk_data.csv': step, commits, shannon, liouville, eff_temp
#   - 'tape_data.csv': step, node, prev, new, entropy, time_stamp
#   - 'sensor_data.csv': step, i, j, heat
# Generates:
#   1. Time-series plots for commits, entropy, Liouville, temperature, and edge heat
#   2. Vortex spiral plot from tape_data
# -----------------------------------------------------------------------------


def load_data():
    df_chunk = pd.read_csv('chunk_data.csv')
    df_tape = pd.read_csv('tape_data.csv')
    df_sensor = pd.read_csv('sensor_data.csv')
    return df_chunk, df_tape, df_sensor


def plot_time_series(df_chunk, df_sensor):
    # Create a 3x2 grid of subplots
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))

    # Commits over time
    ax = axes[0, 0]
    ax.plot(df_chunk['step'], df_chunk['commits'], linewidth=1)
    ax.set_title('Total Commits Over Time')
    ax.set_xlabel('Step'); ax.set_ylabel('Commits'); ax.grid(True)

    # Shannon entropy
    ax = axes[0, 1]
    ax.plot(df_chunk['step'], df_chunk['shannon'], linewidth=1)
    ax.set_title('Shannon Entropy Over Time')
    ax.set_xlabel('Step'); ax.set_ylabel('Entropy'); ax.grid(True)

    # Liouville measure
    ax = axes[1, 0]
    ax.plot(df_chunk['step'], df_chunk['liouville'], linewidth=1)
    ax.set_title('Liouville Measure Over Time')
    ax.set_xlabel('Step'); ax.set_ylabel('Measure'); ax.grid(True)

    # Effective temperature
    ax = axes[1, 1]
    ax.plot(df_chunk['step'], df_chunk['eff_temp'], linewidth=1)
    ax.set_title('Effective Temperature Over Time')
    ax.set_xlabel('Step'); ax.set_ylabel('Temperature (K)'); ax.grid(True)

    # Sensor heat: average heat per step
    df_avg_heat = df_sensor.groupby('step')['heat'].mean().reset_index()
    ax = axes[2, 0]
    ax.plot(df_avg_heat['step'], df_avg_heat['heat'], linewidth=1)
    ax.set_title('Average Edge Heat Over Time')
    ax.set_xlabel('Step'); ax.set_ylabel('Heat'); ax.grid(True)

    # Hide unused subplot
    axes[2, 1].axis('off')

    plt.tight_layout()
    plt.show()


def plot_vortex(df_tape, alpha=1.0):
    # Compute polar to Cartesian for vortex
    theta = 2 * np.pi * (df_tape['entropy'] % 1)
    r = alpha * np.sqrt(df_tape['step'])
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Scatter by node color
    plt.figure(figsize=(6, 6))
    sc = plt.scatter(x, y, c=df_tape['node'], cmap='tab10', s=8)
    plt.title('Quantum Mesh Vortex Spiral')
    plt.xlabel('X'); plt.ylabel('Y')
    plt.axis('equal')
    cbar = plt.colorbar(sc, label='Node ID')
    plt.show()


def plot_transitions(df_tape, alpha=1.0):
    # For each node, connect successive points
    theta = 2 * np.pi * (df_tape['entropy'] % 1)
    r = alpha * np.sqrt(df_tape['step'])
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    df_tape = df_tape.assign(x=x, y=y)

    plt.figure(figsize=(6, 6))
    for node_id, group in df_tape.groupby('node'):
        plt.plot(group['x'], group['y'], linewidth=0.5)
    plt.title('Node Trajectories in Vortex Space')
    plt.xlabel('X'); plt.ylabel('Y')
    plt.axis('equal')
    plt.show()


if __name__ == '__main__':
    df_chunk, df_tape, df_sensor = load_data()
    plot_time_series(df_chunk, df_sensor)
    plot_vortex(df_tape)
    plot_transitions(df_tape)
