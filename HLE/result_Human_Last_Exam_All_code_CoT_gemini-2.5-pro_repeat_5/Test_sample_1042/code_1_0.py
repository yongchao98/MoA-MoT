import pandas as pd
import numpy as np

def solve():
    """
    This function demonstrates that while E(Y^a | A=a, L) is identifiable,
    the full function E(Y^a | A, L) is not because the components where A is not equal to a
    are not identifiable due to unmeasured confounding.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # Define the number of individuals in our simulated population
    n_samples = 20000

    # 1. Simulate the unmeasured confounder U and the measured confounder L
    # U is an unobserved binary factor (e.g., genetic predisposition)
    U = np.random.binomial(1, 0.5, n_samples)
    # L is an observed discrete factor with 3 levels (e.g., age group)
    L = np.random.choice([0, 1, 2], n_samples, p=[0.3, 0.5, 0.2])

    # 2. Simulate the binary treatment A, which is caused by both L and U.
    # P(A=1 | L, U) is a function of L and U, representing self-selection into treatment.
    prob_A_is_1 = 1 / (1 + np.exp(-(-1 + 0.8 * L + 1.5 * U)))
    A = np.random.binomial(1, prob_A_is_1, n_samples)

    # 3. Simulate the outcome Y, which is caused by A, L, and U.
    # The true effect of A on Y is 2.5. U has a strong effect on Y.
    Y = 2.5 * A + 1.5 * L + 3 * U + np.random.normal(0, 2, n_samples)

    # 4. Create the observed dataset. Note that U is not included.
    observed_data = pd.DataFrame({'A': A, 'L': L, 'Y': Y})

    # 5. Let's try to identify E(Y^a | A, L) for the counterfactual a=1.
    # According to our derivation, we can only identify the part where A=1.
    # This is E(Y^{a=1} | A=1, L), which equals E(Y | A=1, L) by consistency.

    print("--- Analysis for Counterfactual a = 1 ---")
    print("The full function E(Y^{a=1} | A, L) is not identifiable.")
    print("However, we can identify the part where A=1.\n")
    
    # Calculate E(Y | A=1, L) from the observed data. This is identifiable.
    identifiable_part = observed_data[observed_data['A'] == 1].groupby('L')['Y'].mean().reset_index()
    identifiable_part.rename(columns={'Y': 'E(Y|A=1,L)'}, inplace=True)

    print("Identifiable Component: E(Y^{a=1} | A=1, L)")
    print("Calculated as the average of observed Y for subjects with A=1, grouped by L:")
    print(identifiable_part.to_string(index=False))
    
    print("\nThe other component, E(Y^{a=1} | A=0, L), is NOT identifiable from the data because")
    print("of the unmeasured confounder U.")

solve()