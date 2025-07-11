import pandas as pd
import numpy as np

def run_simulation():
    """
    Simulates a causal system to demonstrate non-identifiability.
    """
    # 1. Define the Structural Causal Model (SCM) and generate data
    np.random.seed(42)
    n_samples = 1_000_000

    # U is the unmeasured confounder
    U = np.random.binomial(1, 0.5, n_samples)
    # L is the measured confounder
    L = np.random.binomial(1, 0.5, n_samples)

    # A (treatment) is caused by L and U
    # The probability of A=1 depends on L and U
    p_A = 1 / (1 + np.exp(-( -1.5 + 1 * L + 2 * U)))
    A = np.random.binomial(1, p_A, n_samples)

    # Y (outcome) is caused by A, L, and U
    # The true effect of A on Y is 2
    Y = 2 * A + 3 * L + 4 * U + np.random.normal(0, 1, n_samples)

    df = pd.DataFrame({'U': U, 'L': L, 'A': A, 'Y': Y})

    # 2. Define counterfactual outcomes based on the SCM
    # Y^1 is the outcome if everyone received treatment A=1
    # Y^0 is the outcome if everyone received treatment A=0
    df['Y1'] = 2 * 1 + 3 * df['L'] + 4 * df['U'] + np.random.normal(0, 1, n_samples)
    df['Y0'] = 2 * 0 + 3 * df['L'] + 4 * df['U'] + np.random.normal(0, 1, n_samples)

    print("--- Causal Simulation Results ---")
    print("We want to identify E(Y^a | A, L). Let's test for a=1 and L=1.\n")

    # 3. Analyze E(Y^1 | A=1, L=1) - The Identifiable Case
    print("Case 1: Identify E(Y^1 | A=1, L=1)")
    # From the analyst's perspective (using only observed Y, A, L)
    # This relies on the consistency rule: E(Y^1|A=1,L=1) = E(Y|A=1,L=1)
    analyst_view = df[(df['A'] == 1) & (df['L'] == 1)]['Y'].mean()

    # From the oracle's perspective (using the true counterfactual Y^1)
    oracle_view = df[(df['A'] == 1) & (df['L'] == 1)]['Y1'].mean()
    
    print(f"  Analyst's calculation based on observed data E(Y|A=1, L=1): {analyst_view:.4f}")
    print(f"  Oracle's true value E(Y^1|A=1, L=1): {oracle_view:.4f}")
    print("  Conclusion: The values match. This quantity is identifiable.\n")

    # 4. Analyze E(Y^1 | A=0, L=1) - The Non-Identifiable Case
    print("Case 2: Try to identify E(Y^1 | A=0, L=1)")
    # This is the expected outcome under treatment for the group that was NOT treated.
    # The analyst CANNOT compute this directly. We can find the true value with our oracle view.
    oracle_truth = df[(df['A'] == 0) & (df['L'] == 1)]['Y1'].mean()
    print(f"  The true value (from oracle's view) of E(Y^1 | A=0, L=1) is: {oracle_truth:.4f}")
    print("  Can an analyst compute this value from observed data (A, L, Y)?\n")

    # Let's check quantities the analyst COULD compute or might be given.
    # Possibility 1: The analyst naively uses the outcome of the treated group.
    # This is E(Y | A=1, L=1), which we already found is E(Y^1 | A=1, L=1).
    print(f"  Analyst's candidate 1: E(Y|A=1, L=1) = {analyst_view:.4f}. This is not the same.")

    # Possibility 2: The analyst uses the observed outcome of the untreated group.
    # This is E(Y | A=0, L=1), which is E(Y^0 | A=0, L=1) by consistency.
    analyst_y_a0 = df[(df['A'] == 0) & (df['L'] == 1)]['Y'].mean()
    print(f"  Analyst's candidate 2: E(Y|A=0, L=1) = {analyst_y_a0:.4f}. This is not the same.")

    # Possibility 3: The analyst uses the GIVEN identifiable quantity E(Y^1 | L=1).
    # We can calculate this with our oracle view.
    given_quantity = df[df['L'] == 1]['Y1'].mean()
    print(f"  Analyst's candidate 3 (given info): E(Y^1|L=1) = {given_quantity:.4f}. This is not the same.")
    
    print("\n  Conclusion: None of the identifiable quantities match the true value.")
    print("  This is because the subpopulation with A=0 is systematically different")
    print("  from the subpopulation with A=1 due to the unmeasured confounder U.\n")

    print("--- Final Conclusion ---")
    print("Since E(Y^a | A, L) cannot be computed for all values of A, it is not identifiable.")

if __name__ == '__main__':
    run_simulation()