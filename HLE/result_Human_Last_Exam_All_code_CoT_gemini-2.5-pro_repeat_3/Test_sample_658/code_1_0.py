import numpy as np

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value for the minimum degree d_{B,delta} of a polynomial p
    approximating e^{-x} on [0, B] with error < delta.

    The asymptotic formula used is A(B, L) = sqrt(B*L) + L/log(L+e),
    where L = log(1/delta).
    """
    if not (B >= 1 and 0 < delta < 1):
        print("Error: The inputs must satisfy B >= 1 and 0 < delta < 1.")
        return

    L = np.log(1 / delta)

    # Use a safe version of log(L) for the denominator to handle L in (0, inf).
    # log(L + e) is well-defined and positive for L > 0.
    log_L_safe = np.log(L + np.e)
    
    # Calculate the two terms of the asymptotic formula
    term1 = np.sqrt(B * L)
    term2 = L / log_L_safe
    
    # The asymptotic value is the sum of the two terms
    asymptotic_d = term1 + term2
    
    print("This script calculates the asymptotic value for d_{B,delta}.")
    print(f"Given parameters: B = {B}, delta = {delta}")
    print("-" * 20)
    print("Intermediate values:")
    print(f"L = log(1/delta) = {L:.4f}")
    print(f"log(L + e) = {log_L_safe:.4f}")
    print("-" * 20)
    print("Final equation for the asymptotic value A(B, L):")
    print("A(B, L) = sqrt(B * L) + L / log(L + e)")
    print("-" * 20)
    print("Component values from the equation:")
    print(f"Value of 'sqrt(B * L)': {term1:.4f}")
    print(f"Value of 'L / log(L + e)': {term2:.4f}")
    print("-" * 20)
    print(f"The resulting asymptotic value for d_{{B,delta}} is: {asymptotic_d:.4f}")

# You can call this function with your desired values for B and delta.
# For example:
# calculate_asymptotic_degree(100, 0.001)
# calculate_asymptotic_degree(10, 1e-10)

# Let's run an example as requested by the user prompt.
calculate_asymptotic_degree(50, 1e-8)