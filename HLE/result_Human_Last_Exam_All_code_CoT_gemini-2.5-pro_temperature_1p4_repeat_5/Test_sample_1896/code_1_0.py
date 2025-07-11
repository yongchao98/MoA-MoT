import numpy as np
from scipy.integrate import quad

def integrand(x, n):
    """
    The function under the integral for I_n.
    This is the product of n sinc functions.
    """
    # At x=0, the limit of sin(t)/t is 1, so the product is 1.
    if x == 0:
        return 1.0
    
    prod = 1.0
    for k in range(1, n + 1):
        # We use the mathematical definition of sinc(t) = sin(t)/t.
        prod *= np.sin(x / k) / (x / k)
    return prod

def analyze_borwein_integrals(max_n=8):
    """
    Calculates Borwein integrals I_n and related properties to help evaluate the statements.
    """
    pi_half = np.pi / 2
    print("--- Analysis of Borwein Integrals I_n ---")
    print(f"Reference value of π/2 = {pi_half:.15f}\n")
    print(" n | Condition (Σ_2..n 1/k < 1) |    Calculated I_n    | |I_n - π/2|  | P(n) True? ")
    print("---|------------------------------|------------------------|--------------|------------")
    
    # S_n is the sum from k=2 to n of 1/k
    current_sum = 0.0
    for n in range(1, max_n + 1):
        if n >= 2:
            current_sum += 1.0 / n
        
        # Check the condition for I_n = π/2
        condition_holds = current_sum < 1
        if n == 1:
            condition_str = f"S_1=0.0000 (True) "
        else:
            condition_str = f"S_{n}={current_sum:<.4f} ({condition_holds}) "

        # Calculate the integral I_n numerically. np.inf represents infinity.
        # The quad function returns the integral result and an estimated error.
        val, err = quad(integrand, 0, np.inf, args=(n,))
        
        # Check P(n): is I_n = π/2? We use a small tolerance for comparison.
        is_pi_half = abs(val - pi_half) < 1e-9

        # Print the results for this n in a formatted table
        diff = abs(val - pi_half)
        print(f" {n:d} | {condition_str:<28} | {val:.18f} | {diff:1.4e} | {is_pi_half}")

# Run the analysis function
analyze_borwein_integrals()
