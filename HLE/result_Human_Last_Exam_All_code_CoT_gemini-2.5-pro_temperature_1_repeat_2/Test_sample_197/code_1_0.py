import numpy as np

def analyze_correlation_vs_causation():
    """
    Analyzes the relationship between A and D in the system E->A and E->D.
    
    This function simulates data for the structural equation model to demonstrate
    that correlation does not imply causation in the presence of a confounder.
    """
    # 1. Define parameters for the simulation
    num_samples = 100000  # Number of data points
    np.random.seed(42)    # for reproducibility

    # Causal coefficients representing the strength of the causal links
    # E -> A relationship equation: A = (coef_e_to_a * E) + noise_a
    # E -> D relationship equation: D = (coef_e_to_d * E) + noise_d
    coef_e_to_a = 0.8
    coef_e_to_d = -0.7

    # 2. Generate data based on the structural equations
    # E is an exogenous variable (a common cause)
    e_data = np.random.randn(num_samples)
    
    # A and D are caused by E, plus some random noise.
    # There is NO direct causal link between A and D in these equations.
    noise_a = np.random.randn(num_samples)
    a_data = coef_e_to_a * e_data + noise_a

    noise_d = np.random.randn(num_samples)
    d_data = coef_e_to_d * e_data + noise_d

    # 3. Calculate the correlation between A and D
    correlation_matrix = np.corrcoef(a_data, d_data)
    correlation_ad = correlation_matrix[0, 1]

    # 4. Explain the result and print the final answer
    print("Structural Equation Model Simulation:")
    print(f"  - Equation for A: A = {coef_e_to_a}*E + noise_A")
    print(f"  - Equation for D: D = {coef_e_to_d}*E + noise_D")
    print("\nIn this model, A and D are not causally related to each other.")
    print("However, they share a common cause, E.")
    print("\nLet's calculate the correlation from the simulated data:")
    print(f"  - The calculated correlation between A and D is: {correlation_ad:.4f}")
    
    print("\nConclusion:")
    print("The correlation is strong and non-zero. This is because the common cause 'E' creates a spurious association between A and D.")
    print("Therefore, in this system, correlation does not imply causation.")
    print("\nDoes correlation imply causation in this system?")
    print("No")

analyze_correlation_vs_causation()