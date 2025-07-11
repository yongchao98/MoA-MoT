import numpy as np
import pandas as pd

def analyze_causation():
    """
    Simulates the given structural equation model to determine if correlation
    implies causation between variables A and D.
    """
    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # Number of data points to generate for the observational study
    n_samples = 10000

    # Step 1: Define the structural equations and generate data
    # The model is E->A->B->C<-D<-E
    # We will represent this with linear equations and random noise.

    # E is an exogenous variable (a common cause for A and D)
    # Equation: E = N(0, 2^2) where N is a normal distribution
    e_std = 2
    E = np.random.normal(loc=0, scale=e_std, size=n_samples)

    # A is caused by E
    # Equation: A = 1.5 * E + N(0, 1^2)
    a_from_e = 1.5
    A = a_from_e * E + np.random.normal(loc=0, scale=1, size=n_samples)

    # D is also caused by E
    # Equation: D = -1.0 * E + N(0, 1^2)
    d_from_e = -1.0
    D = d_from_e * E + np.random.normal(loc=0, scale=1, size=n_samples)

    # B is caused by A
    # Equation: B = 0.8 * A + N(0, 1^2)
    b_from_a = 0.8
    B = b_from_a * A + np.random.normal(loc=0, scale=1, size=n_samples)

    # C is caused by B and D (making C a collider)
    # Equation: C = 1.2 * B + 0.9 * D + N(0, 1^2)
    c_from_b = 1.2
    c_from_d = 0.9
    C = c_from_b * B + c_from_d * D + np.random.normal(loc=0, scale=1, size=n_samples)

    # Create a DataFrame to hold the simulated data
    df = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E})

    # Step 2: Calculate the correlation between A and D
    correlation_A_D = df['A'].corr(df['D'])

    # Step 3: Explain the results and provide the final answer
    print("--- Causal System Simulation ---")
    print("Model: E->A->B->C<-D<-E")
    print(f"\nSimulating {n_samples} data points based on the following equations:")
    print(f"E = Normal(mean=0, std={e_std})")
    print(f"A = {a_from_e} * E + Normal(mean=0, std=1)")
    print(f"D = {d_from_e} * E + Normal(mean=0, std=1)")
    print(f"B = {b_from_a} * A + Normal(mean=0, std=1)")
    print(f"C = {c_from_b} * B + {c_from_d} * D + Normal(mean=0, std=1)")

    print("\n--- Analysis of Correlation vs. Causation ---")
    print(f"The calculated correlation between A and D is: {correlation_A_D:.4f}")
    print("\nA strong correlation between A and D is observed in the data.")
    print("However, in the model we defined, there is no direct or indirect causal path from A to D or D to A.")
    print("This correlation exists because both A and D share a common cause, E (A <- E -> D).")
    print("This is a classic example of confounding, where an external variable creates a spurious correlation.")
    print("\nTherefore, in this system, the correlation between A and D does not imply causation.")

    print("\n<<<No>>>")

if __name__ == '__main__':
    analyze_causation()