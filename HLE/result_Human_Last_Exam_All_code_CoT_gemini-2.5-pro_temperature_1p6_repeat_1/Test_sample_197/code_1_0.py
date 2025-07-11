import numpy as np

def simulate_and_analyze_sem():
    """
    Simulates the Structural Equation Model (SEM) E->A->B->C<-D<-E
    and checks the correlation between A and D.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # Number of data points to simulate
    n_samples = 10000

    # Define the causal relationships with random coefficients.
    # The key is that E is a common cause of A and D.
    # E -> A
    # E -> D
    c_ea = 0.8
    c_ed = 0.7
    
    # We can model the other relationships as well, though they are not
    # needed to show the A-D correlation.
    # A -> B
    # B -> C
    # D -> C
    c_ab = 0.6
    c_bc = 0.5
    c_dc = 0.4

    # Generate the base variable E and noise terms from a standard normal distribution
    E = np.random.randn(n_samples)
    noise_a = np.random.randn(n_samples)
    noise_d = np.random.randn(n_samples)
    noise_b = np.random.randn(n_samples)
    noise_c = np.random.randn(n_samples)
    
    # Generate the variables based on the structural equations
    # The model equations are:
    # A = c_ea * E + noise_a
    # D = c_ed * E + noise_d
    A = c_ea * E + noise_a
    D = c_ed * E + noise_d
    B = c_ab * A + noise_b
    C = c_bc * B + c_dc * D + noise_c

    # Calculate the correlation between A and D
    # The correlation matrix is returned, so we take the value at (0, 1)
    correlation_ad = np.corrcoef(A, D)[0, 1]

    # --- Output the analysis ---
    print("Analysis of the Structural Equation Model: E->A->B->C<-D<-E\n")
    print(f"The system is simulated with {n_samples} samples.")
    print("In this model, there is no direct or indirect causal path between A and D.")
    print("However, both A and D share a common cause: E.")
    print(f"The equation for A is: A = {c_ea} * E + noise")
    print(f"The equation for D is: D = {c_ed} * E + noise")
    print("\nThis common cause (confounder) induces a correlation between A and D.")
    
    print("\n--- Simulation Result ---")
    print(f"The calculated correlation between A and D is: {correlation_ad:.4f}")

    print("\n--- Conclusion ---")
    print("Since the correlation is caused by the confounding variable E and not by a causal link between A and D,")
    print("the correlation does not imply causation in this system.")


# Run the simulation and print the results
simulate_and_analyze_sem()