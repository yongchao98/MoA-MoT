import argparse

def calculate_critical_correlation(N, mu, theta):
    """
    Calculates the critical amount of correlation required to balance
    potentiation and depression in the described neural network model.

    Args:
        N (int): The number of neurons in each input layer (N_v = N_s = N).
        mu (float): The average firing rate of input neurons.
        theta (float): The heterosynaptic offset constant (plasticity threshold).
    """

    # For a Poisson process, the variance of the rate is equal to the mean rate.
    sigma_sq = mu

    # Derived formula for the critical covariance C = Cov(v_k, s_k)
    # C = 2 * N * mu * (theta - mu) - sigma^2
    critical_covariance = 2 * N * mu * (theta - mu) - sigma_sq

    # The correlation coefficient rho is the covariance divided by the product of standard deviations.
    # Since Var(v) = Var(s) = sigma_sq, rho = C / sigma_sq
    correlation_coefficient = critical_covariance / sigma_sq if sigma_sq != 0 else float('inf')


    print("Calculating the critical covariance 'C' based on the derived formula:")
    print(f"C = 2*N*μ*(θ - μ) - σ²")
    print("\nGiven the parameters:")
    print(f"  N (neurons per layer) = {N}")
    print(f"  μ (mean rate)         = {mu}")
    print(f"  θ (threshold)         = {theta}")
    print(f"  σ² (variance, assumed μ) = {sigma_sq}")

    print("\nThe calculation is:")
    print(f"C = 2 * {N} * {mu} * ({theta} - {mu}) - {sigma_sq}")
    print(f"C = {2*N*mu} * ({theta-mu}) - {sigma_sq}")
    print(f"C = {2*N*mu*(theta-mu)} - {sigma_sq}")
    
    print("\n----------------------------------------------------")
    print(f"Resulting critical covariance: C = {critical_covariance}")
    print(f"Corresponding correlation coefficient: ρ = {correlation_coefficient:.4f}")
    print("----------------------------------------------------")

    if not -1.0 <= correlation_coefficient <= 1.0:
        print("\nWarning: The calculated correlation coefficient is outside the physical range of [-1, 1].")
        print("This suggests that for the given parameters (N, μ, θ), no amount of correlation can")
        print("balance potentiation and depression in this model.")
    else:
        print("\nThe calculated correlation coefficient is within the physical range [-1, 1].")


if __name__ == '__main__':
    # You can run this script from the command line with arguments,
    # or modify the default values below and run it directly.
    parser = argparse.ArgumentParser(description="Calculate critical correlation in a neural network.")
    parser.add_argument('--N', type=int, default=100, help='Number of neurons per layer.')
    parser.add_argument('--mu', type=float, default=10.0, help='Average firing rate (e.g., in Hz).')
    parser.add_argument('--theta', type=float, default=10.04, help='Plasticity threshold (in the same units as mu).')
    
    args = parser.parse_args()
    
    calculate_critical_correlation(args.N, args.mu, args.theta)

    # Final answer format for the derived formula
    # Let C be the critical amount of correlation (covariance).
    # The result is the formula expressed in terms of the given parameters.
    # C = 2*N*μ*(θ - μ) - σ²
    # Assuming the input is a Poisson process, variance σ² = μ.
    # So, C = 2*N*μ*(θ - μ) - μ
    final_formula = "2*N*μ*(θ - μ) - μ"
    print(f"\n<<<C = {final_formula}>>>")
