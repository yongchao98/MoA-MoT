import numpy as np

def simulate_and_correlate():
    """
    Simulates data based on the SEM to demonstrate confounding.
    """
    # for reproducibility of the random data
    np.random.seed(42)

    # Number of data points in the observational study
    n_samples = 10000

    # Define the causal coefficients for the model's equations
    # These represent the strength of the causal links
    coef_E_to_A = 0.8
    coef_E_to_D = -0.7 # Using a negative coefficient for variety

    # Simulate the data based on the structural equations
    # E is an exogenous variable (a source), modeled as a standard normal distribution
    E = np.random.normal(loc=0, scale=1, size=n_samples)

    # A and D are caused by E, plus some independent random noise
    noise_A = np.random.normal(loc=0, scale=1, size=n_samples)
    noise_D = np.random.normal(loc=0, scale=1, size=n_samples)

    # The final structural equations for A and D.
    # Note that there is no term for A in the equation for D, or vice-versa.
    # This means there is no direct causation between them in our model.
    A = coef_E_to_A * E + noise_A
    D = coef_E_to_D * E + noise_D

    # We don't need to simulate B and C for this problem as they are not
    # involved in the relationship between A, D, and their common cause E.

    # Calculate the Pearson correlation between the simulated A and D variables
    correlation_matrix = np.corrcoef(A, D)
    correlation_AD = correlation_matrix[0, 1]

    # Output the results
    print("A simulation of the system: E->A and E->D")
    print("------------------------------------------")
    print(f"The structural equation for A is: A = {coef_E_to_A} * E + noise")
    print(f"The structural equation for D is: D = {coef_E_to_D} * E + noise")
    print("\nAs shown in the equations, there is no direct causal link between A and D.")
    print("However, they share a common cause, E, which creates a spurious correlation.\n")
    print(f"Based on a simulation with {n_samples} data points:")
    # Here, we output the "final equation" which is the correlation calculation result
    print(f"The calculated correlation between A and D is: {correlation_AD:.4f}")

if __name__ == '__main__':
    simulate_and_correlate()