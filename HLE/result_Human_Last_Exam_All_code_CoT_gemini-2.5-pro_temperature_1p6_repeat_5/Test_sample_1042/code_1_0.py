import pandas as pd
import numpy as np

def solve_causal_identification():
    """
    Demonstrates the identification of E(Y^a | A, L) given E(Y^a | L).

    The script follows these steps:
    1.  Sets up a Structural Causal Model (SCM) with an unmeasured confounder U.
    2.  Generates a large dataset including observed variables (A, L, Y)
        and unobserved counterfactuals (Y^0, Y^1).
    3.  Specifies a target: Identify E(Y^1 | A, L) for the stratum L=1.
    4.  Simulates an "oracle" providing the value of E(Y^1 | L=1), which is a premise
        of the problem. In a real-world scenario, this would come from another source
        or method.
    5.  Uses the identification formula to compute E(Y^1 | A=0, L=1) from
        observational data and the oracle value.
    6.  Compares the identified value with the ground truth from the simulation.
    """
    np.random.seed(42)
    n_samples = 200000

    # 1. Define SCM and generate data
    # U is the unmeasured confounder
    U = np.random.binomial(1, 0.5, n_samples)
    # L is a measured confounder that is also affected by U
    L = np.random.binomial(1, 0.2 + 0.5 * U)
    # A is the treatment, affected by L and U
    A = np.random.binomial(1, 0.1 + 0.2 * L + 0.4 * U)
    # Y is the outcome
    Y = 1.5 * A + 0.8 * L + 2.0 * U + np.random.normal(0, 1, n_samples)

    # Generate counterfactual outcomes
    Y0 = 1.5 * 0 + 0.8 * L + 2.0 * U + np.random.normal(0, 1, n_samples)
    Y1 = 1.5 * 1 + 0.8 * L + 2.0 * U + np.random.normal(0, 1, n_samples)

    df = pd.DataFrame({'U': U, 'L': L, 'A': A, 'Y': Y, 'Y0': Y0, 'Y1': Y1})

    # --- Analysis for stratum L=1 ---
    print("--- Analysis for Stratum L=1 and Counterfactual a=1 ---")
    df_l1 = df[df['L'] == 1]

    # 2. Get the "oracle" value for E(Y^a | L)
    # This is given by the problem's premise. We calculate it from our
    # full dataset (including Y^1) to simulate the oracle.
    # E[Y^a|L] is identified.
    oracle_E_Y1_L1 = df_l1['Y1'].mean()

    # 3. Use the identification strategy on the "observed" data
    # We want to identify E(Y^1 | A, L=1).
    # This requires identifying E(Y^1 | A=1, L=1) and E(Y^1 | A=0, L=1).

    # Case A=a (i.e., A=1, a=1): Identify via consistency
    # E(Y^1 | A=1, L=1) = E(Y | A=1, L=1)
    observed_E_Y1_A1_L1 = df_l1[df_l1['A'] == 1]['Y'].mean()

    # Get probabilities P(A|L=1) from observed data
    P_A1_L1 = len(df_l1[df_l1['A'] == 1]) / len(df_l1)
    P_A0_L1 = 1 - P_A1_L1

    # Case A!=a (i.e., A=0, a=1): Identify using the formula
    print("\nCalculating E(Y^1 | A=0, L=1) using the identification formula:")
    print("Formula: [E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1)] / P(A=0|L=1)\n")
    print(f"Value for E(Y^1|L=1) (from oracle): {oracle_E_Y1_L1:.4f}")
    print(f"Value for E(Y^1|A=1,L=1) (identified from data): {observed_E_Y1_A1_L1:.4f}")
    print(f"Value for P(A=1|L=1) (from data): {P_A1_L1:.4f}")
    print(f"Value for P(A=0|L=1) (from data): {P_A0_L1:.4f}")

    # The equation with numbers plugged in:
    numerator = oracle_E_Y1_L1 - observed_E_Y1_A1_L1 * P_A1_L1
    identified_E_Y1_A0_L1 = numerator / P_A0_L1
    
    print(f"\nFinal Equation:")
    print(f"E(Y^1 | A=0, L=1) = ({oracle_E_Y1_L1:.4f} - {observed_E_Y1_A1_L1:.4f} * {P_A1_L1:.4f}) / {P_A0_L1:.4f} = {identified_E_Y1_A0_L1:.4f}")


    # 4. Compare with the ground truth
    # We can calculate the true values directly from our simulation because we have Y^1
    true_E_Y1_A0_L1 = df_l1[df_l1['A'] == 0]['Y1'].mean()
    true_E_Y1_A1_L1 = df_l1[df_l1['A'] == 1]['Y1'].mean() # same as observed_E_Y1_A1_L1 by consistency

    print("\n--- Verification ---")
    print("We compare the identified values with the ground truth from the simulation.")
    print("\nIdentified value for E(Y^1 | A=1, L=1): {:.4f}".format(observed_E_Y1_A1_L1))
    print("Ground truth for E(Y^1 | A=1, L=1):   {:.4f}".format(true_E_Y1_A1_L1))
    print("\nIdentified value for E(Y^1 | A=0, L=1): {:.4f}".format(identified_E_Y1_A0_L1))
    print("Ground truth for E(Y^1 | A=0, L=1):   {:.4f}".format(true_E_Y1_A0_L1))
    
    # Check the condition E(Y^a | A,L) != E(Y^a | L)
    print("\nFinally, we check the problem's condition that E(Y^a|A,L) != E(Y^a|L):")
    print(f"E(Y^1|L=1) = {oracle_E_Y1_L1:.4f}")
    print(f"E(Y^1|A=0,L=1) = {true_E_Y1_A0_L1:.4f}")
    print(f"E(Y^1|A=1,L=1) = {true_E_Y1_A1_L1:.4f}")
    print("The values are different, so the condition holds.")

solve_causal_identification()