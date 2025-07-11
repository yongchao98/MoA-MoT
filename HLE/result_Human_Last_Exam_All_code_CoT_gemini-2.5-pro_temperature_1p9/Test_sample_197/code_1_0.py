import numpy as np
import pandas as pd

def simulate_data_and_show_correlation():
    """
    Simulates data based on the structural equation model E->A->B->C<-D<-E
    and calculates the correlation between A and D.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # Number of data points
    n_samples = 2000

    # 1. Generate data based on the causal graph
    # E is an exogenous variable (has no causes in the model)
    E = np.random.randn(n_samples)

    # A and D are caused by E
    # A = 0.7*E + noise
    A = 0.7 * E + np.random.randn(n_samples)
    # D = -0.6*E + noise
    D = -0.6 * E + np.random.randn(n_samples)

    # B is caused by A
    # B = 0.5*A + noise
    B = 0.5 * A + np.random.randn(n_samples)

    # C is a "collider," caused by B and D
    # C = 0.4*B + 0.8*D + noise
    C = 0.4 * B + 0.8 * D + np.random.randn(n_samples)

    # 2. Create a pandas DataFrame for easier analysis
    data = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E})

    # 3. Calculate the correlation matrix
    correlation_matrix = data.corr()

    # 4. Extract and print the correlation between A and D
    correlation_A_D = correlation_matrix.loc['A', 'D']
    
    print("This simulation creates data where A and D are not causally related,")
    print("but they share a common cause E.")
    print("The resulting correlation between A and D is calculated below.")
    print("-" * 50)
    print(f"Calculated Correlation(A, D): {correlation_A_D:.4f}")
    print("-" * 50)
    print("As you can see, the correlation is not zero, which is caused by the confounding variable E.")
    print("This demonstrates that correlation does not imply causation in this system.")


simulate_data_and_show_correlation()
