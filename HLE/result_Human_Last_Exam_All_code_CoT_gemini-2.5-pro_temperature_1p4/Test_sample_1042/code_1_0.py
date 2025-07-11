import numpy as np
import pandas as pd

def sigmoid(x):
    """Sigmoid function for converting log-odds to probability."""
    return 1 / (1 + np.exp(-x))

def run_simulation_and_analysis():
    """
    Simulates a causal system with an unmeasured confounder and demonstrates
    the non-identifiability of E(Y^a | A!=a, L).
    """
    # Set seed for reproducibility
    np.random.seed(42)
    
    # 1. Define a Causal Model and Generate "True" Data
    # A large population size
    n_samples = 1_000_000

    # U: Unmeasured confounder (binary)
    U = np.random.binomial(1, 0.5, n_samples)
    
    # L: Measured confounder (binary)
    L = np.random.binomial(1, 0.4, n_samples)
    
    # A: Treatment, depends on L and U
    # logit(P(A=1)) = -1 + 2*L + 2*U
    prob_A_is_1 = sigmoid(-1 + 2 * L + 2 * U)
    A = np.random.binomial(1, prob_A_is_1, n_samples)
    
    # Y^a: Potential outcomes, depend on L and U. Y^1 is outcome if A=1, Y^0 if A=0.
    # The true treatment effect is 5.
    Y0 = 10 + 3 * L + 6 * U + np.random.normal(0, 1, n_samples)
    Y1 = 15 + 3 * L + 6 * U + np.random.normal(0, 1, n_samples)
    
    # Y: Observed outcome, determined by the consistency rule (Y = Y^A)
    Y = Y1 * A + Y0 * (1 - A)
    
    # Create a full "God's-eye-view" dataframe
    full_df = pd.DataFrame({'U': U, 'L': L, 'A': A, 'Y0': Y0, 'Y1': Y1, 'Y': Y})

    # 2. Calculate the "True" Quantities of Interest
    # We want to know E(Y^a | A, L). Let's calculate these for L=0.
    # E(Y^a | A=a', L=l)
    
    # --- When a = a' ---
    # E[Y^1 | A=1, L=0]: The expected outcome under treatment for those treated, with L=0
    true_Y1_given_A1_L0 = full_df[(full_df['A'] == 1) & (full_df['L'] == 0)]['Y1'].mean()
    # E[Y^0 | A=0, L=0]: The expected outcome under control for those in control, with L=0
    true_Y0_given_A0_L0 = full_df[(full_df['A'] == 0) & (full_df['L'] == 0)]['Y0'].mean()
    
    # --- When a != a' ---
    # E[Y^1 | A=0, L=0]: The expected outcome under treatment for those in control, with L=0
    true_Y1_given_A0_L0 = full_df[(full_df['A'] == 0) & (full_df['L'] == 0)]['Y1'].mean()
    # E[Y^0 | A=1, L=0]: The expected outcome under control for those treated, with L=0
    true_Y0_given_A1_L0 = full_df[(full_df['A'] == 1) & (full_df['L'] == 0)]['Y0'].mean()

    # 3. Emulate the Analyst with Observed Data
    # The analyst only sees (A, L, Y).
    observed_df = full_df[['A', 'L', 'Y']]
    
    # Calculate what is possible from observed data
    # The analyst can calculate E[Y | A, L]
    obs_Y_given_A1_L0 = observed_df[(observed_df['A'] == 1) & (observed_df['L'] == 0)]['Y'].mean()
    obs_Y_given_A0_L0 = observed_df[(observed_df['A'] == 0) & (observed_df['L'] == 0)]['Y'].mean()

    # 4. Print and Compare the Results
    
    print("This script demonstrates the non-identifiability of E(Y^a | A, L) due to an unmeasured confounder U.")
    print("We will focus on the subpopulation where the measured confounder L=0.\n")

    print("--- Analysis for E(Y^1 | A, L=0) ---")
    print("Goal: Find the average outcome if everyone had treatment 1, for different groups.")
    
    print(f"\nCase 1: (A=1) What is E(Y^1 | A=1, L=0)? (Counterfactual treatment matches observed treatment)")
    print(f"    - True value from full data: E[Y^1 | A=1, L=0] = {true_Y1_given_A1_L0:.3f}")
    print(f"    - Analyst's calculation from observed data: E[Y | A=1, L=0] = {obs_Y_given_A1_L0:.3f}")
    print("    - Conclusion: The values match. This quantity is identifiable because Y=Y^1 for this group.")

    print(f"\nCase 2: (A=0) What is E(Y^1 | A=0, L=0)? (Counterfactual treatment differs from observed treatment)")
    print(f"    - True value from full data: E[Y^1 | A=0, L=0] = {true_Y1_given_A0_L0:.3f}")
    print(f"    - Analyst's calculation for this group is for Y, which is Y^0 not Y^1: E[Y | A=0, L=0] = {obs_Y_given_A0_L0:.3f}")
    print("    - Conclusion: The analyst has no way to calculate E[Y^1 | A=0, L=0]. Their observed data for this group, E[Y|A=0,L=0], relates to Y^0, not Y^1.")
    print(f"    - The true value ({true_Y1_given_A0_L0:.3f}) is different from what the analyst observes for this group ({obs_Y_given_A0_L0:.3f}) and also different from the identifiable E(Y^1|A=1,L=0) ({true_Y1_given_A1_L0:.3f}).")
    print("\nBecause E(Y^a|A,L) is not identifiable for A != a, the quantity is not generally identifiable.")


if __name__ == '__main__':
    run_simulation_and_analysis()