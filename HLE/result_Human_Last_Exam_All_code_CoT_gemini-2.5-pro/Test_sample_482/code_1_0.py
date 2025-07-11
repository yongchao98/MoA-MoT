import math

def calculate_critical_correlation(mu, N, tau_r):
    """
    Calculates the critical correlation C required to balance potentiation and depression.

    Args:
        mu (float): The average rate of activation for input neurons (in Hz).
        N (int): The number of neurons in each input population.
        tau_r (float): The leaky integrator time constant for output neurons (in seconds).
    
    Returns:
        float: The critical covariance C (in Hz^2).
    """
    
    # The variance of the input rate is sigma^2 = mu / tau_r
    sigma_squared = mu / tau_r
    
    # The critical covariance C is sigma^2 / N
    critical_covariance = sigma_squared / N
    
    return critical_covariance

def main():
    """
    Main function to define parameters and print the result.
    """
    # Placeholder parameters for the model
    # mu: average firing rate of input neurons (e.g., 10 Hz)
    mu = 10.0
    # N: number of neurons in each input layer (e.g., 100)
    N = 100
    # tau_r: membrane time constant of output neurons (e.g., 20ms = 0.02s)
    tau_r = 0.02

    print("Calculating the critical amount of correlation (C)...")
    print("The formula is: C = mu / (N * tau_r)")
    print("-" * 30)
    
    # Calculate the result
    c_critical = calculate_critical_correlation(mu, N, tau_r)
    
    # Print the equation with the numbers plugged in
    print(f"Given parameters:")
    print(f"  mu (average rate) = {mu} Hz")
    print(f"  N (neurons per layer) = {N}")
    print(f"  tau_r (time constant) = {tau_r} s")
    print("-" * 30)

    print("Final Equation:")
    print(f"C = {mu} / ({N} * {tau_r})")
    
    # Print the final calculated value
    print(f"C = {c_critical} Hz^2")
    
    # The final answer in the requested format
    print("\n<<<" + str(c_critical) + ">>>")

if __name__ == "__main__":
    main()
