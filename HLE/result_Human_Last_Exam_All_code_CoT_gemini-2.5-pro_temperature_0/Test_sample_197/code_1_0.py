import numpy as np
import pandas as pd
from scipy.stats import pearsonr

def simulate_sem_and_test_correlation():
    """
    Simulates the Structural Equation Model (SEM) E->A->B->C<-D<-E
    and tests the correlation between A and D.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # Number of data points in the observational study
    n_samples = 1000

    # 1. Define the structural equations based on the model
    # We assume linear relationships for simplicity.
    # E is an exogenous variable (a root cause), so it's just random noise.
    e = np.random.randn(n_samples)

    # A and D are caused by E.
    # A = 1.5*E + noise
    a = 1.5 * e + np.random.randn(n_samples)
    # D = -1.2*E + noise
    d = -1.2 * e + np.random.randn(n_samples)

    # B is caused by A.
    # B = 0.8*A + noise
    b = 0.8 * a + np.random.randn(n_samples)

    # C is caused by B and D (C is a collider).
    # C = 0.6*B + 0.9*D + noise
    c = 0.6 * b + 0.9 * d + np.random.randn(n_samples)

    # 2. Analyze the relationship between A and D
    # The data was generated with no direct causal link between A and D.
    # A is caused by E. D is also caused by E.
    # This common cause (confounder) E should induce a correlation.

    # Calculate the Pearson correlation between A and D
    correlation, p_value = pearsonr(a, d)

    # 3. Print the results
    print("Simulating the system: E -> A -> B -> C <- D <- E")
    print("-" * 50)
    print(f"Number of samples: {n_samples}")
    print("In our simulation, the structural equations are:")
    print("E = N(0, 1)")
    print("A = 1.5*E + N(0, 1)")
    print("D = -1.2*E + N(0, 1)")
    print("\nNote: There is no direct causal term for 'A' in the equation for 'D', or vice-versa.")
    print("-" * 50)
    print(f"The calculated correlation between A and D is: {correlation:.4f}")
    print(f"The p-value is: {p_value}")
    print("\nConclusion: A strong correlation is found between A and D.")
    print("However, this correlation is due to the common cause 'E', not a causal relationship between A and D.")
    print("Therefore, in this system, correlation does not imply causation.")

# Run the simulation
simulate_sem_and_test_correlation()

# Final answer as a single word
print("\nDoes correlation imply causation in this system?")
print("<<<No>>>")