import numpy as np

def calculate_critical_correlation():
    """
    This function calculates the critical correlation between two input populations (v and s)
    required to balance potentiation and depression in a feedforward network.
    """

    # --- Model Parameters (Example Values) ---
    # Number of neurons in each input layer (assuming N_v = N_s = N)
    N = 100
    # Average firing rate of input neurons (e.g., in Hz)
    mu = 5.0
    # Heterosynaptic depression threshold (e.g., in Hz)
    theta = 5.015
    # Variance of the input neuron firing rates (e.g., in Hz^2)
    sigma_in_sq = 10.0

    print("--- Calculating the Critical Amount of Correlation ---")
    print("\nWe are given the following parameters:")
    print(f"  Number of inputs (N_v=N_s): N = {N}")
    print(f"  Average input rate: mu = {mu}")
    print(f"  Depression threshold: theta = {theta}")
    print(f"  Input rate variance: sigma_in^2 = {sigma_in_sq}")

    # For a physically possible correlation to exist (i.e., correlation coefficient <= 1),
    # the parameters must satisfy the condition: N*mu*(theta-mu) <= sigma_in^2.
    # This condition ensures that the required covariance is not greater than the variance.
    depression_term = N * mu * (theta - mu)

    print(f"\nFirst, we check the solvability condition: N*mu*(theta-mu) <= sigma_in^2")
    print(f"  {N} * {mu} * ({theta} - {mu}) <= {sigma_in_sq}")
    print(f"  {depression_term:.4f} <= {sigma_in_sq}")

    if depression_term > sigma_in_sq:
        print("\nResult: Condition NOT met.")
        print("For these parameters, depression is too strong to be balanced by potentiation,")
        print("even with maximum positive correlation between the inputs.")
        return

    print("\nResult: Condition met. A physical solution exists.")

    # Calculate the critical covariance (C_vs)
    # The formula is derived from a stability analysis of the weight dynamics:
    # C_vs = 2 * N * mu * (theta - mu) - sigma_in^2
    critical_C_vs = 2 * depression_term - sigma_in_sq

    # Calculate the critical correlation coefficient (rho)
    critical_rho = critical_C_vs / sigma_in_sq

    print("\nThe critical covariance (C_vs) required to balance potentiation and depression is given by the formula:")
    print("  C_vs = 2 * N * mu * (theta - mu) - sigma_in^2")
    print("\nPlugging in the numbers:")
    print(f"  C_vs = 2 * {N} * {mu} * ({theta} - {mu}) - {sigma_in_sq}")
    print(f"  C_vs = 2 * {N * mu} * ({theta - mu:.3f}) - {sigma_in_sq}")
    print(f"  C_vs = {2 * depression_term} - {sigma_in_sq}")
    print(f"  C_vs = {critical_C_vs:.4f}")

    print("\nThe corresponding dimensionless correlation coefficient (rho) is:")
    print("  rho = C_vs / sigma_in^2")
    print(f"  rho = {critical_C_vs:.4f} / {sigma_in_sq}")
    print(f"  rho = {critical_rho:.4f}")

if __name__ == '__main__':
    calculate_critical_correlation()
