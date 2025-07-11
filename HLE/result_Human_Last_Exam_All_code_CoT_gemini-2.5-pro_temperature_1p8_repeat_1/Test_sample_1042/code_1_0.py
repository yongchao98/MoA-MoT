import pandas as pd
import numpy as np

def demonstrate_identification_issue(n_samples=200000, seed=42):
    """
    Simulates a causal system with an unmeasured confounder to test identification.

    Causal Structure:
    - U (unmeasured) and L (measured) are confounders.
    - U, L -> A (treatment)
    - U, L, A -> Y (outcome)
    """
    np.random.seed(seed)
    
    # U is an unmeasured confounder (e.g., genetic predisposition)
    U = np.random.binomial(1, 0.5, n_samples)
    # L is a measured confounder (e.g., age group)
    L = np.random.binomial(1, 0.4, n_samples)
    
    # A is treatment assignment, dependent on both L and U
    # The effect of U is stronger to create a significant confounding bias
    A_prob = 1 / (1 + np.exp(-(0.7 * L + 1.8 * U - 1.2)))
    A = np.random.binomial(1, A_prob, n_samples)
    
    # Y^a are the potential outcomes. Both are affected by U and L.
    # Let Y^1 be the outcome if treated (A=1) and Y^0 if not treated (A=0).
    # The treatment A has a causal effect of +5 on Y.
    Y0 = 2.0 * L + 3.0 * U + np.random.normal(0, 1, n_samples)
    Y1 = 5.0 + 2.0 * L + 3.0 * U + np.random.normal(0, 1, n_samples)
    
    # The observed outcome Y is the potential outcome corresponding to the actual treatment received.
    Y = A * Y1 + (1 - A) * Y0
    
    # This is the "God's view" or "true" dataset, including unobservables
    full_data = pd.DataFrame({'U': U, 'L': L, 'A': A, 'Y': Y, 'Y0': Y0, 'Y1': Y1})
    
    # Let's focus our analysis on the stratum where the measured confounder L=1
    L_val = 1
    
    # =========================================================================
    # Step 1: Using the "true" data, find the ground truth values we want to identify.
    # The quantity is E(Y^a | A, L). We check for a=1 and L=1.
    
    # Case A=a (i.e., A=1): The true value of E(Y^1 | A=1, L=1)
    true_val_A_is_a = full_data[(full_data['A'] == 1) & (full_data['L'] == L_val)]['Y1'].mean()
    
    # Case A!=a (i.e., A=0): The true value of E(Y^1 | A=0, L=1)
    true_val_A_is_not_a = full_data[(full_data['A'] == 0) & (full_data['L'] == L_val)]['Y1'].mean()
    # =========================================================================

    # =========================================================================
    # Step 2: Now, use only the "observed" data a researcher would have: (A, L, Y).
    observed_data = full_data[['A', 'L', 'Y']]

    # Attempt to identify E(Y^1 | A=1, L=1) from observed data.
    # By consistency, E(Y^1 | A=1, L=1) = E(Y | A=1, L=1). This is identifiable.
    identified_val = observed_data[(observed_data['A'] == 1) & (observed_data['L'] == L_val)]['Y'].mean()

    # Attempt to identify E(Y^1 | A=0, L=1).
    # A researcher CANNOT compute this from observed data, as they never see Y1 for this group.
    # Any estimate would require an untestable assumption. For example, assuming no unmeasured confounding
    # would mean E(Y^1 | A=0, L=1) = E(Y^1 | A=1, L=1), which is false here.
    # =========================================================================
    
    print("--- Goal: Can we identify E(Y^a | A, L)? ---")
    print("Let's test for a=1 and L=1.")
    print("\nCASE 1: A=a. We want to identify E(Y^1 | A=1, L=1).")
    print("This is the expected outcome under treatment for the group that was actually treated.")
    print(f"The TRUE value is E[Y^1 | A=1, L=1] = {true_val_A_is_a:.4f}")
    print(f"The IDENTIFIED value using observed data is E[Y | A=1, L=1] = {identified_val:.4f}")
    print("Result: Success! Identification is possible when A=a.")
    
    print("\nCASE 2: A!=a. We want to identify E(Y^1 | A=0, L=1).")
    print("This is the expected outcome under treatment for the group that was NOT treated.")
    print(f"The TRUE value is E[Y^1 | A=0, L=1] = {true_val_A_is_not_a:.4f}")
    print("Result: Failure! We cannot calculate this from observed data (A, L, Y).")
    print(f"If we naively used the value from the treated group ({identified_val:.4f}), we would be wrong.")
    print(f"The bias is {abs(true_val_A_is_not_a - identified_val):.4f}, caused by the unmeasured confounder U.")
    
    print("\n--- Final Conclusion ---")
    print("Because E(Y^a | A, L) cannot be computed for all its arguments (specifically when A!=a),")
    print("the entire quantity is considered NOT identifiable.")

demonstrate_identification_issue()