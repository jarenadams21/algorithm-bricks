import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.errors import EmptyDataError

# —— 1) TRY TO LOAD output.csv FOR R PLOTS ——
df_out = None
if os.path.exists('output.csv'):
    try:
        df_out = pd.read_csv('output.csv')
        # drop any empty columns
        df_out = df_out.dropna(axis=1, how='all')
        # strip whitespace
        df_out.columns = df_out.columns.str.strip()
        # rename if necessary
        if 'R Bell' in df_out.columns:
            df_out = df_out.rename(columns={'R Bell': 'R'})
    except EmptyDataError:
        print("output.csv is empty — skipping entanglement‐R plots")
else:
    print("output.csv not found — skipping entanglement‐R plots")

# —— 2) LOAD summary.csv (must exist) ——
df_sum = pd.read_csv('summary.csv')
df_sum.columns = df_sum.columns.str.strip()

# —— 3) LOAD attack.csv ——
df_atk = pd.read_csv('attack.csv')
df_atk.columns = df_atk.columns.str.strip()

# Ensure numeric types
for col in ['step','P_mean','S_mean_inside','S_mean_outside']:
    df_sum[col] = pd.to_numeric(df_sum[col], errors='coerce')
for col in ['step','px','py','R','S','inside']:
    if col in df_atk.columns:
        df_atk[col] = pd.to_numeric(df_atk[col], errors='coerce')

# —— PLOTS —— #

# A) entanglement reliability over time (if df_out loaded)
if df_out is not None and 'R' in df_out.columns:
    plt.figure()
    plt.plot(df_out['step'], df_out['R'], label='Entanglement R')
    plt.xlabel('Step')
    plt.ylabel('Reliability')
    plt.title('Entanglement Reliability Over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Phase‐space: R vs P_mean
    if {'P01','P12','P20'}.issubset(df_out.columns):
        df_out['P_mean'] = df_out[['P01','P12','P20']].mean(axis=1)
        plt.figure()
        plt.scatter(df_out['R'], df_out['P_mean'])
        plt.xlabel('Entanglement R')
        plt.ylabel('Mean Route Probability')
        plt.title('Entanglement Phase Space')
        plt.grid(True)
        plt.show()

else:
    print("Skipping entanglement phase‐space plots.")

# B) Transport vs Attack Success from summary.csv
plt.figure()
plt.plot(df_sum['step'], df_sum['P_mean'],      label='Mean Route Probability')
plt.plot(df_sum['step'], df_sum['S_mean_inside'], label='Mean Attack Success (Inside)')
plt.plot(df_sum['step'], df_sum['S_mean_outside'],label='Mean Attack Success (Outside)')
plt.xlabel('Step')
plt.ylabel('Probability')
plt.title('Transport vs Attack Success Over Time')
plt.legend()
plt.grid(True)
plt.show()

# C) Classify inner vs outer attack points
if 'inside' in df_atk.columns:
    inner_mask = df_atk['inside'].astype(bool)
else:
    # fallback: consider S close to (1-R) as inner
    inner_mask = np.isclose(df_atk['S'], 1 - df_atk['R'], atol=1e-8)

# spatial scatter
plt.figure()
plt.scatter(df_atk.loc[inner_mask, 'px'],
            df_atk.loc[inner_mask, 'py'], alpha=0.7)
plt.xlabel('px'); plt.ylabel('py')
plt.title('Attack Points Inside Mesh')
plt.grid(True)
plt.show()

plt.figure()
plt.scatter(df_atk.loc[~inner_mask, 'px'],
            df_atk.loc[~inner_mask, 'py'], alpha=0.7)
plt.xlabel('px'); plt.ylabel('py')
plt.title('Attack Points Outside Mesh')
plt.grid(True)
plt.show()

# boxplot of S
plt.figure()
plt.boxplot([df_atk.loc[inner_mask, 'S'],
             df_atk.loc[~inner_mask, 'S']],
            labels=['Inner','Outer'])
plt.ylabel('Success Probability S')
plt.title('Attack Success Distribution')
plt.grid(True)
plt.show()
