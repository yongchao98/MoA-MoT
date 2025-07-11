import math

def calculate_kappa(mu, tau_v, phi, rho, A):
    """
    Calculates the parameter kappa based on the derived formula.

    Args:
        mu (float): The mean firing rate of synapses.
        tau_v (float): The time constant for the presynaptic accumulator.
        phi (float): The scaling constant for the presynaptic accumulator.
        rho (float): The offset constant in the Hebbian equation.
        A (float): The power/variance of the input signal x_k.

    Returns:
        float: The value of kappa.
    """
    numerator = -mu * tau_v * (phi * mu + rho)
    denominator = phi * A
    
    if denominator == 0:
        return float('inf') if numerator > 0 else float('-inf') if numerator < 0 else float('nan')
        
    kappa = numerator / denominator
    return kappa

def main():
    """
    Main function to define parameters, calculate kappa, and print the result.
    """
    # Example physiological parameters
    mu = 10.0      # Firing rate in Hz
    tau_v = 0.020  # Presynaptic time constant in seconds (20 ms)
    phi = 1.0      # Presynaptic scaling constant
    rho = -15.0    # Offset constant (LTD-inducing)
    
    # For a Poisson process, the power spectral density is equal to the rate.
    A = mu         # Power of the input signal
    
    # Calculate kappa
    kappa = calculate_kappa(mu, tau_v, phi, rho, A)
    
    # Print the definition and the result with numbers
    print("The definition of kappa is: kappa = - (mu * tau_v * (phi * mu + rho)) / (phi * A)")
    print("\nUsing the following parameter values:")
    print(f"mu = {mu}")
    print(f"tau_v = {tau_v}")
    print(f"phi = {phi}")
    print(f"rho = {rho}")
    print(f"A = {A}")

    print("\nWe calculate kappa = - ({} * {} * ({} * {} + {})) / ({} * {})".format(
        mu, tau_v, phi, mu, rho, phi, A
    ))
    
    # The final format is specified by the user prompt
    print("\n<<<{}>>>".format(kappa))

if __name__ == "__main__":
    main()