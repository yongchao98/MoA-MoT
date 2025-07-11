import numpy as np
import pandas as pd

def check_correlation_vs_causation():
    """
    Analyzes and simulates the SEM E->A->B->C<-D<-E to determine if
    the correlation between A and D implies causation.
    """

    # --- Step 1: Theoretical Analysis ---
    print("Analyzing the causal model: E->A->B->C<-D<-E")
    print("-------------------------------------------------")
    print("In this model, the variable E is a direct cause of both A and D.")
    print("This makes E a 'common cause' or 'confounder' for A and D.")
    print("A common cause can create a statistical correlation between its effects (A and D) without any direct or indirect causal link between them.")
    print("Therefore, an observed correlation between A and D is expected due to the confounding effect of E.")
    print("\n--- Step 2: Computational Simulation ---")
    print("To demonstrate this, we will simulate data based on the model.")

    # --- Model Parameters ---
    n_samples = 10000
    # Coefficients for the causal links
    coef_e_a = 0.8
    coef_e_d = -0.7 # Using a negative coefficient to show correlation still occurs
    coef_a_b = 0.6
    coef_b_c = 0.5
    coef_d_c = 0.4
    # Standard deviation for the random noise
    noise_std = 0.3

    # --- Data Generation ---
    # E is an exogenous variable (it has no causes in this model)
    E = np.random.randn(n_samples)
    
    # A and D are caused by E, plus some random noise
    A = coef_e_a * E + np.random.normal(0, noise_std, n_samples)
    D = coef_e_d * E + np.random.normal(0, noise_std, n_samples)

    # B is caused by A, C is caused by B and D
    # These are included for model completeness but are not needed for the A-D correlation calculation
    B = coef_a_b * A + np.random.normal(0, noise_std, n_samples)
    C = coef_b_c * B + coef_d_c * D + np.random.normal(0, noise_std, n_samples)
    
    # Put data into a pandas DataFrame for easy correlation calculation
    df = pd.DataFrame({'A': A, 'B': B, 'C': C, 'D': D, 'E': E})

    # --- Step 3: Calculate and Report Correlation ---
    print("\nThe structural equations used for A and D are:")
    print(f"A = {coef_e_a} * E + noise")
    print(f"D = {coef_e_d} * E + noise")

    # Calculate the correlation between A and D
    correlation_ad = df['A'].corr(df['D'])

    print(f"\nAfter simulating {n_samples} observations, the calculated correlation between A and D is: {correlation_ad:.4f}")
    
    # --- Step 4: Final Conclusion ---
    print("\n--- Conclusion ---")
    print("The simulation shows a strong correlation between A and D.")
    print("However, we know from the model's structure that this is not due to a causal relationship between them.")
    print("It is a spurious correlation caused by the confounder, E.")
    print("\nTherefore, does correlation imply causation in this system?")
    print("No.")

if __name__ == '__main__':
    check_correlation_vs_causation()