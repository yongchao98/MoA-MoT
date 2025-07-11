import numpy as np
import pandas as pd

def simulate_and_analyze():
    """
    Simulates the causal system E->A->B->C<-D<-E and checks the correlation between A and D.
    """
    # Define the number of data points for the simulation
    n_samples = 20000

    # Define the causal coefficients for the structural equations
    # These numbers represent the strength of the causal influence.
    coef_e_a = 0.8
    coef_e_d = 0.7
    coef_a_b = 0.6
    coef_b_c = 0.5
    coef_d_c = 0.4

    # 1. Simulate E, the common cause. It's an exogenous variable (no causes in this model).
    # E ~ Normal(0, 1)
    e = np.random.randn(n_samples)

    # 2. Simulate A and D. Both are caused by E.
    # The final equations are:
    # A = 0.8 * E + noise
    # D = 0.7 * E + noise
    noise_a = np.random.randn(n_samples) * 0.2
    a = coef_e_a * e + noise_a

    noise_d = np.random.randn(n_samples) * 0.2
    d = coef_e_d * e + noise_d
    
    # 3. Simulate B and C based on their causes
    # B = 0.6 * A + noise
    noise_b = np.random.randn(n_samples) * 0.2
    b = coef_a_b * a + noise_b

    # C = 0.5 * B + 0.4 * D + noise
    noise_c = np.random.randn(n_samples) * 0.2
    c = coef_b_c * b + coef_d_c * d + noise_c

    # 4. Create a DataFrame to hold the simulated data
    df = pd.DataFrame({'A': a, 'B': b, 'C': c, 'D': d, 'E': e})

    # 5. Calculate the correlation between A and D
    correlation_ad = df['A'].corr(df['D'])

    # 6. Print the results and conclusion
    print("--- Causal System Simulation ---")
    print("Model: E->A->B->C<-D<-E")
    print("\nSimulated Structural Equations (with coefficients):")
    print(f"A = {coef_e_a}*E + Noise")
    print(f"D = {coef_e_d}*E + Noise")
    print(f"B = {coef_a_b}*A + Noise")
    print(f"C = {coef_b_c}*B + {coef_d_c}*D + Noise")
    
    print("\n--- Analysis ---")
    print("As shown in the equations above, 'A' is not a cause of 'D', and 'D' is not a cause of 'A'.")
    print(f"\nHowever, after simulating {n_samples} data points, the calculated correlation is:")
    print(f"Correlation(A, D) = {correlation_ad:.4f}")

    print("\nThis strong, positive correlation is generated because 'E' is a common cause of both 'A' and 'D'.")
    print("This is a classic example of confounding, where a non-causal correlation is observed.")

if __name__ == '__main__':
    simulate_and_analyze()

<<<No>>>