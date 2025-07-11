import pandas as pd
import numpy as np

def simulate_data_and_show_correlation():
    """
    Simulates data from the structural equation model E->A->B->C<-D<-E
    and calculates the correlation between A and D.
    """
    # Set a seed for reproducibility
    np.random.seed(42)
    
    # Number of data points to generate
    n_samples = 10000
    
    # Define the causal relationships with linear equations and random noise.
    # The coefficients represent the strength of the causal effect.
    # E is an exogenous variable (it has no specified cause in the model)
    E = np.random.randn(n_samples)
    
    # A is caused by E
    # A = 1.0 * E + noise
    A = 1.0 * E + np.random.randn(n_samples) * 0.5
    
    # D is also caused by E
    # D = 0.8 * E + noise
    D = 0.8 * E + np.random.randn(n_samples) * 0.5
    
    # B is caused by A
    # B = 0.7 * A + noise
    B = 0.7 * A + np.random.randn(n_samples) * 0.5
    
    # C is caused by B and D (C is a collider)
    # C = 0.6 * B + 0.9 * D + noise
    C = 0.6 * B + 0.9 * D + np.random.randn(n_samples) * 0.5
    
    # Create a pandas DataFrame to hold the simulated data
    data = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E})
    
    # Calculate the correlation between A and D
    correlation_A_D = data['A'].corr(data['D'])
    
    print("Simulated Structural Equation Model: E->A->B->C<-D<-E")
    print("In this model, A and D share a common cause E, but do not cause each other.")
    print(f"Based on {n_samples} simulated data points:")
    print(f"The calculated correlation between A and D is: {correlation_A_D:.4f}")
    print("\nThis non-zero correlation is induced by the common cause E.")
    print("Therefore, in this system, correlation does not imply causation.")
    print("\nFinal Answer:")
    # The final answer to the question "Does correlation imply causation in this system?"
    # is a single word.
    final_answer = "No"
    print(final_answer)

# Run the simulation
simulate_data_and_show_correlation()