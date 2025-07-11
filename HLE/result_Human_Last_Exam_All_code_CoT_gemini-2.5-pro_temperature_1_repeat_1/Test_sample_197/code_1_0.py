import pandas as pd
import numpy as np

def simulate_data(n_samples=10000):
    """
    Simulates data for the structural equation model:
    E -> A -> B -> C <- D <- E
    """
    # Set random seed for reproducibility
    np.random.seed(42)

    # 1. E is an exogenous variable (a common cause for A and D)
    # We can model it as a standard normal distribution.
    e = np.random.randn(n_samples)

    # 2. A is caused by E. D is caused by E.
    # A = 0.8 * E + noise
    # D = 0.7 * E + noise
    # The coefficients (0.8, 0.7) represent the strength of the causal effect.
    a = 0.8 * e + np.random.randn(n_samples) * 0.5
    d = 0.7 * e + np.random.randn(n_samples) * 0.5

    # 3. B is caused by A.
    # B = 0.6 * A + noise
    b = 0.6 * a + np.random.randn(n_samples) * 0.5

    # 4. C is caused by B and D (C is a collider)
    # C = 0.5 * B + 0.4 * D + noise
    c = 0.5 * b + 0.4 * d + np.random.randn(n_samples) * 0.5

    # Create a pandas DataFrame to hold the data
    df = pd.DataFrame({'A': a, 'B': b, 'C': c, 'D': d, 'E': e})
    return df

def analyze_correlation(df):
    """
    Calculates and prints the correlation between A and D.
    """
    # Calculate the correlation matrix
    correlation_matrix = df.corr()

    # Get the specific correlation between A and D
    correlation_ad = correlation_matrix.loc['A', 'D']

    print("Simulated Structural Equation Model: E->A->B->C<-D<-E")
    print("\nIn this model, A and D share a common cause, E, but have no direct causal relationship.")
    print("\nCalculating the correlation between A and D from the simulated data:")
    print(f"The correlation between A and D is: {correlation_ad:.4f}")
    print("\nThis non-zero correlation is induced by the confounding variable E.")
    print("Therefore, the observed correlation does not imply causation between A and D.")


# --- Main execution ---
# Simulate the data based on the model
simulated_df = simulate_data()

# Analyze the correlation
analyze_correlation(simulated_df)
