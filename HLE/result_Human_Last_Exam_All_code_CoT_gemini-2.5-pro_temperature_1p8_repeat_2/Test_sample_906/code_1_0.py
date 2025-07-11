import math

def calculate_steady_state_pi0(rho):
    """
    Calculates the steady-state probability pi_0 for a birth-death process.
    The formula derived is pi_0 = e^(-rho).

    Args:
        rho (float): The traffic intensity, defined as lambda / mu.

    Returns:
        float: The value of pi_0.
    """
    return math.exp(-rho)

def main():
    """
    Main function to demonstrate the calculation.
    """
    # Example values for lambda and mu
    # Let's assume lambda = 2 and mu = 1, so rho = 2.0
    lambda_rate = 2.0
    mu_rate = 1.0
    rho = lambda_rate / mu_rate

    # The final derived formula for pi_0 is e**(-rho)
    # The numbers in this final equation are 0 (in pi_0) and the exponent, which is -rho.
    # In python terms, that's math.exp(-rho).
    
    print("The derived formula for the steady-state probability pi_0 is:")
    print("pi_0 = e**(-rho)")
    
    # Calculate pi_0 for the example rho
    pi_0 = calculate_steady_state_pi0(rho)

    print(f"\nFor an example where rho = {rho}:")
    # Output the final equation with the specific numbers
    print(f"pi_0 = e**(-{rho})")
    
    print(f"\nThe calculated value of pi_0 is: {pi_0}")

if __name__ == "__main__":
    main()
