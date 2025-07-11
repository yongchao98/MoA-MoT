import numpy as np
import pandas as pd

def analyze_causal_model():
    """
    Simulates data based on the SEM E->A->B->C<-D<-E to check if correlation
    between A and D implies causation.
    """
    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # Number of data points in the observational study
    n_samples = 2000

    # 1. Define the structural equations based on the model.
    # We will use simple linear relationships plus some random noise.
    # The equations are:
    # E = noise_E
    # A = 1.0*E + noise_A
    # D = 0.8*E + noise_D
    # B = 0.7*A + noise_B
    # C = 0.6*B + 0.9*D + noise_C
    
    print("Simulating data based on the following structural equations:")
    print("E = N(0, 1)")
    print("A = 1.0 * E + N(0, 1)")
    print("D = 0.8 * E + N(0, 1)")
    # The other variables are generated for completeness of the model, but not
    # essential for analyzing the A-D relationship.
    print("B = 0.7 * A + N(0, 1)")
    print("C = 0.6 * B + 0.9 * D + N(0, 1)")
    print("-" * 30)

    # 2. Generate data for each variable
    # E is an exogenous variable (no causes in this graph)
    E = np.random.randn(n_samples)
    
    # A is caused by E
    A = 1.0 * E + np.random.randn(n_samples)
    
    # D is caused by E
    D = 0.8 * E + np.random.randn(n_samples)

    # B is caused by A
    B = 0.7 * A + np.random.randn(n_samples)
    
    # C is a "collider", caused by both B and D
    C = 0.6 * B + 0.9 * D + np.random.randn(n_samples)
    
    # 3. Create a pandas DataFrame to hold the data
    df = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E})
    
    # 4. Calculate the correlation between A and D
    # In our simulation, there is no direct causal arrow A->D or D->A.
    # However, they share a common cause E. This should create a correlation.
    correlation_A_D = df['A'].corr(df['D'])
    
    # 5. Print the result
    # We expect a non-zero correlation, showing that correlation exists.
    # But we know from our model setup that there is no causation between them.
    print(f"The simulated correlation between A and D is: {correlation_A_D:.4f}")
    print("\nThis correlation is not due to a causal link between A and D.")
    print("It is a spurious correlation caused by the common variable 'E' influencing both 'A' and 'D'.")
    print("Therefore, in this system, correlation does not imply causation.")

if __name__ == '__main__':
    analyze_causal_model()
