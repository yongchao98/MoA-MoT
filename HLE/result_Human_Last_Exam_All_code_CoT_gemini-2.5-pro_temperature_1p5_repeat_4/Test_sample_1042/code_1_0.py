import pandas as pd
import numpy as np

def run_simulation_and_identification():
    """
    Simulates data with an unmeasured confounder and demonstrates
    the identification of E(Y^a | A, L).
    """
    # 1. Define a Structural Causal Model (SCM) and generate data
    np.random.seed(42)
    n_samples = 200000
    # L is a measured confounder, U is an unmeasured confounder
    L = np.random.binomial(1, 0.5, n_samples)
    U = np.random.binomial(1, 0.4, n_samples)
    
    # A is caused by L and U
    A_prob = 1 / (1 + np.exp(-(L + 2*U - 1.5))) 
    A = np.random.binomial(1, A_prob, n_samples)
    
    # Y is caused by A, L, and U
    noise = np.random.normal(0, 1, n_samples)
    # Define counterfactual outcomes Y^0 and Y^1
    Y0 = 1.5 * L + 2 * U + noise
    Y1 = 2 + 1.5 * L + 2 * U + noise
    
    # Observed outcome Y depends on the actual treatment A
    Y = A * Y1 + (1 - A) * Y0

    # Create a DataFrame with all variables
    # In a real scenario, U, Y0, Y1 are unobserved ('God's view')
    full_data = pd.DataFrame({'L': L, 'U': U, 'A': A, 'Y': Y, 'Y0': Y0, 'Y1': Y1})
    # This is the dataset an analyst would observe
    observed_data = full_data[['A', 'L', 'Y']]

    # We want to identify E(Y^a | A, L) for a=1
    a = 1

    print("--- Ground Truth Calculation (using unobserved U and Y^a) ---")
    # Group by L and A to compute the true E(Y^a|A,L)
    ground_truth = full_data.groupby(['L', 'A'])[f'Y{a}'].mean().reset_index()
    ground_truth = ground_truth.rename(columns={f'Y{a}': 'True E(Y^a|A,L)'})
    print(ground_truth)
    print("\n--------------------------------------------------------------\n")
    
    print("--- Identification from Observed Data ---")
    
    # According to the problem, E(Y^a | L) is identified. For this simulation,
    # we "identify" it by computing it from the full data, then provide it
    # to the main calculation as a known quantity.
    # E_Ya_given_L is a function of (a,L) computable from P(Y,A,L).
    E_Ya_given_L = full_data.groupby('L')[f'Y{a}'].mean()
    print("Assumed identified quantity E(Y^a|L) for a=1:")
    print(E_Ya_given_L)
    print("\n")
    
    # Let's perform the identification for each level of L (0 and 1)
    for l_val in [0, 1]:
        print(f"--- Calculations for L = {l_val} ---")
        
        # --- Part 1: Identify E(Y^a | A=a, L) ---
        print(f"1. Identifying E(Y^a | A=a, L={l_val}) for a={a}:")
        
        # We use E(Y^a | A=a, L) = E(Y | A=a, L)
        obs_subset = observed_data[(observed_data['A'] == a) & (observed_data['L'] == l_val)]
        identified_E_Ya_given_A_eq_a = obs_subset['Y'].mean()
        
        print(f"   E(Y^a | A={a}, L={l_val}) is identified by E(Y | A={a}, L={l_val})")
        print(f"   Value = {identified_E_Ya_given_A_eq_a:.4f}\n")
        
        # --- Part 2: Identify E(Y^a | A=1-a, L) ---
        print(f"2. Identifying E(Y^a | A=1-a, L={l_val}) for a={a}:")

        # Observational probabilities (propensity scores conditional on L)
        p_A_eq_a = len(observed_data[(observed_data['A'] == a) & (observed_data['L'] == l_val)]) / len(observed_data[observed_data['L'] == l_val])
        p_A_neq_a = 1 - p_A_eq_a

        # The identified value of E(Y^a|L) for this l_val
        val_E_Ya_given_L = E_Ya_given_L.loc[l_val]

        # Use the formula: E(Y^a | A=1-a, L) = [E(Y^a|L) - E(Y|A=a,L) * P(A=a|L)] / P(A=1-a|L)
        numerator = val_E_Ya_given_L - (identified_E_Ya_given_A_eq_a * p_A_eq_a)
        identified_E_Ya_given_A_neq_a = numerator / p_A_neq_a

        print(f"   Using formula: [E(Y^a|L) - E(Y|A=a,L) * P(A=a|L)] / P(A=1-a|L)")
        print(f"   Values: [{val_E_Ya_given_L:.4f} - {identified_E_Ya_given_A_eq_a:.4f} * {p_A_eq_a:.4f}] / {p_A_neq_a:.4f}")
        print(f"   Numerator = {numerator:.4f}")
        print(f"   Value = {identified_E_Ya_given_A_neq_a:.4f}\n")
    
        print(f"--- Summary for L = {l_val} ---")
        true_val_A_eq_a = ground_truth.loc[(ground_truth['L'] == l_val) & (ground_truth['A'] == a), 'True E(Y^a|A,L)'].iloc[0]
        true_val_A_neq_a = ground_truth.loc[(ground_truth['L'] == l_val) & (ground_truth['A'] == 1-a), 'True E(Y^a|A,L)'].iloc[0]
        print(f"Identified E(Y^{a}|A={a}, L={l_val}): {identified_E_Ya_given_A_eq_a:.4f} \t vs. True Value: {true_val_A_eq_a:.4f}")
        print(f"Identified E(Y^{a}|A={1-a}, L={l_val}): {identified_E_Ya_given_A_neq_a:.4f} \t vs. True Value: {true_val_A_neq_a:.4f}")
        print("-" * 30)

if __name__ == '__main__':
    run_simulation_and_identification()