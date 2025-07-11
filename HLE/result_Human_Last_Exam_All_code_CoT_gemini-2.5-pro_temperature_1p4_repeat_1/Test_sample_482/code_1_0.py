import numpy as np

def calculate_critical_covariance(N, mu, theta):
    """
    Calculates the critical input covariance required to balance potentiation and depression.

    The derivation starts from the weight dynamics dW/dt ~ r(x - theta), where x is the
    input (v or s). For the total weight onto an output neuron r_i to be stable, the
    time-averaged change must be zero.
    
    Condition for stability: E[r_i * (sum(v_k) + sum(s_k) - (N_v + N_s)*theta)] = 0
    Assuming N_v = N_s = N, this becomes:
    Cov(r_i, sum(v_k + s_k)) = (2*N*theta - E[sum(v_k + s_k)]) * E[r_i]
    
    Using the approximations that r_i is mainly driven by v_i and s_i, and that
    inputs v_j, v_k (j!=k) and v_i, s_k (i!=k) are uncorrelated, we can solve for
    the critical covariance between corresponding inputs, C_vs = Cov(v_i, s_i).
    
    The final derived formula is:
    C_vs = mu * (2 * N * (theta - mu) - 1)

    Args:
        N (int): The number of neurons in each input layer (N_v = N_s = N).
        mu (float): The average rate of activation for inputs (v and s).
        theta (float): The heterosynaptic depression constant.

    Returns:
        float: The critical covariance C_vs.
    """

    # --- Define example parameters ---
    # We choose plausible values where the depression threshold is slightly
    # higher than the mean activation rate, creating a net depressive force
    # that needs to be balanced by correlation.
    # N: Number of neurons per layer
    # mu: Average input firing rate (e.g., in Hz)
    # theta: Firing rate threshold for depression (e.g., in Hz)

    # --- Calculation ---
    # This is the derived expression for the critical covariance Cov(v_i, s_i)
    critical_covariance = mu * (2 * N * (theta - mu) - 1)
    
    # --- Outputting the results as requested ---
    print("To find the critical amount of correlation, we derived a formula for the covariance between")
    print("corresponding neurons in the two input populations (Cov(v_i, s_i)) that balances potentiation and depression.")
    print("\nThe derived formula is: C_vs = mu * (2 * N * (theta - mu) - 1)")
    
    print("\nUsing the following example parameter values:")
    print(f"  N (neurons per layer) = {N}")
    print(f"  mu (average input rate) = {mu}")
    print(f"  theta (depression threshold) = {theta}")
    
    print("\nWe can substitute these values into the equation:")
    # Print each part of the equation with numbers
    part1 = "C_vs = " + str(mu)
    part2 = " * (2 * " + str(N)
    part3 = " * (" + str(theta) + " - " + str(mu) + ") - 1)"
    print(part1 + part2 + part3)

    # Print intermediate calculation steps
    term_in_parentheses = 2 * N * (theta - mu) - 1
    print(f"C_vs = {mu} * ({2 * N * (theta - mu)} - 1)")
    print(f"C_vs = {mu} * ({term_in_parentheses})")
    
    print("\nThe critical amount of correlation (covariance) is:")
    print(f"{critical_covariance}")
    
    return critical_covariance

# --- Main execution block ---
if __name__ == "__main__":
    # Example Parameters
    N_neurons = 100
    mu_rate = 0.5
    theta_threshold = 0.6
    
    # Calculate and print the result
    final_covariance = calculate_critical_covariance(N_neurons, mu_rate, theta_threshold)
    # The final answer is wrapped according to the format requirement.
    # This captures the numerical result from the example calculation.
    print(f"\n<<<{final_covariance}>>>")
