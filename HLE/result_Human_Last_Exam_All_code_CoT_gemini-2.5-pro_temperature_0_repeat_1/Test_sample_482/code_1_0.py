import numpy as np

def calculate_critical_correlation():
    """
    Calculates the critical amount of correlation (covariance) required 
    to balance potentiation and depression in the network.
    """
    # --- Parameters ---
    # Number of neurons in the input layers
    N_v = 50
    N_s = 50
    
    # Average rate of activation for input neurons (in Hz)
    mu = 5.0
    
    # Heterosynaptic offset constant (in Hz)
    theta = 5.015

    # --- Calculation ---
    # The formula derived from the stability analysis is:
    # c = -μ + (N_v + N_s) * μ * (θ - μ)
    
    # Calculate the critical covariance 'c'
    c = -mu + (N_v + N_s) * mu * (theta - mu)

    # --- Output ---
    # Print the equation with the numerical values
    print("The equation for the critical covariance 'c' is:")
    print(f"c = -μ + (N_v + N_s) * μ * (θ - μ)")
    print("\nPlugging in the values:")
    print(f"c = -{mu} + ({N_v} + {N_s}) * {mu} * ({theta} - {mu})")
    
    # Print the final result
    print(f"\nThe calculated critical covariance is: {c}")

    # The correlation coefficient rho can also be calculated as c/mu
    # A valid correlation requires |rho| <= 1.
    rho = c / mu
    print(f"This corresponds to a correlation coefficient (rho = c/mu) of: {rho:.3f}")
    if abs(rho) > 1:
        print("\nWarning: The required correlation coefficient is outside the valid range of [-1, 1].")
        print("This means that for the given parameters, the network cannot be stabilized by this correlation alone.")

if __name__ == "__main__":
    calculate_critical_correlation()
