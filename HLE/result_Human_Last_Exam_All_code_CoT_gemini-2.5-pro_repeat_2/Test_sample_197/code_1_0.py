import numpy as np

def simulate_and_analyze():
    """
    Simulates the structural equation model to check if correlation implies causation
    for variables A and D.
    """
    # Define the number of data points for the simulation
    n_samples = 10000

    # Define the causal coefficients for the model E->A and E->D
    # These represent the strength of the causal links.
    e_to_a_coeff = 0.8
    e_to_d_coeff = -0.7 # Using a negative coefficient to show correlation can be negative too

    # 1. Generate data based on the structural equations.
    # In the model E->A<-E->D, E is a common cause for A and D.
    # We represent this with the following equations:
    # A = (coefficient * E) + random_noise
    # D = (coefficient * E) + random_noise

    # Generate the common cause 'E' from a standard normal distribution
    E = np.random.normal(loc=0, scale=1, size=n_samples)

    # Generate 'A' as a function of 'E' plus some random noise
    noise_A = np.random.normal(loc=0, scale=0.5, size=n_samples)
    A = e_to_a_coeff * E + noise_A

    # Generate 'D' as a function of 'E' plus some random noise
    noise_D = np.random.normal(loc=0, scale=0.5, size=n_samples)
    D = e_to_d_coeff * E + noise_D

    # 2. Calculate the correlation between A and D
    # np.corrcoef returns a 2x2 matrix. The correlation between A and D is at [0, 1]
    correlation_ad = np.corrcoef(A, D)[0, 1]

    # 3. Print the results and the explanation
    print("This script simulates the causal structure where A and D share a common cause, E.")
    print("The specific equations used in this simulation are:")
    print(f"A = {e_to_a_coeff} * E + noise")
    print(f"D = {e_to_d_coeff} * E + noise")
    print("-" * 30)
    print(f"Based on {n_samples} generated data points:")
    print(f"The calculated correlation between A and D is: {correlation_ad:.4f}")
    print("\nConclusion from simulation:")
    print("A significant correlation is found between A and D.")
    print("This correlation arises because both variables are caused by E.")
    print("In the model, there is no causal arrow from A to D or from D to A.")
    print("Thus, the correlation does not imply causation in this system.")

if __name__ == '__main__':
    simulate_and_analyze()
