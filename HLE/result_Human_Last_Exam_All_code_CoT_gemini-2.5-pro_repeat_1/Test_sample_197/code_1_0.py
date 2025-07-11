import numpy as np

def simulate_causal_system():
    """
    Simulates the causal system E->A->B->C<-D<-E to check the correlation between A and D.
    """
    # Set a seed for reproducibility of the random numbers
    np.random.seed(42)

    # Number of data points to simulate
    num_samples = 5000

    # Define the causal coefficients for the linear equations
    # These represent the strength of the causal links
    # E -> A
    c_ea = 0.8
    # E -> D
    c_ed = -0.7
    # A -> B
    c_ab = 0.6
    # B -> C
    c_bc = 0.5
    # D -> C
    c_dc = 0.9

    # 1. Generate data for the exogenous variable 'E' (the common cause)
    # Let's assume E is drawn from a standard normal distribution
    e = np.random.randn(num_samples)

    # 2. Generate data for A and D, both caused by E
    # A = c_ea * E + noise_A
    noise_a = np.random.randn(num_samples) * 0.5
    a = c_ea * e + noise_a

    # D = c_ed * E + noise_D
    noise_d = np.random.randn(num_samples) * 0.5
    d = c_ed * e + noise_d
    
    # 3. Generate data for the rest of the system (for completeness)
    # B = c_ab * A + noise_B
    noise_b = np.random.randn(num_samples) * 0.4
    b = c_ab * a + noise_b
    
    # C = c_bc * B + c_dc * D + noise_C
    noise_c = np.random.randn(num_samples) * 0.3
    c = c_bc * b + c_dc * d + noise_c

    # 4. Calculate the correlation between A and D
    correlation_matrix = np.corrcoef(a, d)
    correlation_ad = correlation_matrix[0, 1]

    # Print the parameters of the simulated equations
    print("Simulated Structural Equation Model: E->A->B->C<-D<-E\n")
    print("Causal Equations used in the simulation:")
    print(f"A = {c_ea}*E + noise_A")
    print(f"D = {c_ed}*E + noise_D")
    print(f"B = {c_ab}*A + noise_B")
    print(f"C = {c_bc}*B + {c_dc}*D + noise_C")
    print("-" * 40)

    # Print the result
    print(f"The calculated correlation between A and D is: {correlation_ad:.4f}")
    print("\nConclusion:")
    print("Despite the strong correlation, there is no direct or indirect causal path between A and D in our model.")
    print("The correlation arises because 'E' is a common cause to both A and D.")
    print("This is a classic case where correlation does not imply causation.")

if __name__ == '__main__':
    simulate_causal_system()
