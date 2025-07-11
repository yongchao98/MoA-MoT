import numpy as np
from scipy.special import factorial
from math import comb

def solve_and_print():
    """
    This function implements the plan to solve the problem.
    It calculates the Taylor coefficients, finds the smallest n for which f(n) > 10,
    computes the infinity norm of W_n, and prints the final result.
    """
    MAX_DEGREE = 50  # A reasonable upper bound for the degree n

    # Step 1: Calculate Taylor coefficients for g(x) = (2/pi)K(x)exp(x)
    # The coefficients for (2/pi)K(x) are k_i = (binom(2i, i) / 4^i)^2
    # The coefficients for exp(x) are e_i = 1/i!
    # The final coefficients c_i are the Cauchy product of k_i and e_i.

    k_coeffs = np.zeros(MAX_DEGREE + 1, dtype=float)
    for i in range(MAX_DEGREE + 1):
        if i == 0:
            k_coeffs[i] = 1.0
        else:
            k_coeffs[i] = (comb(2 * i, i) / (4**i))**2

    e_coeffs = 1.0 / factorial(np.arange(MAX_DEGREE + 1), exact=False)

    g_coeffs = np.zeros(MAX_DEGREE + 1, dtype=float)
    for i in range(MAX_DEGREE + 1):
        g_coeffs[i] = np.dot(k_coeffs[:i + 1], e_coeffs[i::-1])

    n_final = -1
    final_roots = None
    final_f_n = 0.0

    # Step 2: Loop to find the smallest n where f(n) > 10
    for n in range(1, MAX_DEGREE + 1):
        # P_n(x) has degree n, coefficients are c_0, ..., c_n
        # np.roots wants coefficients from highest power to lowest
        poly_coeffs = g_coeffs[:n + 1][::-1]

        roots = np.roots(poly_coeffs)

        f_n = np.sum(np.abs(roots)**3)

        if f_n > 10:
            n_final = n
            final_roots = roots
            final_f_n = f_n
            break

    if n_final == -1:
        print("Failed to find the required n within the specified degree limit.")
        return

    # Step 3: Calculate ||W_n||_infinity
    # Handle potential repeated roots for the norm calculation.
    # We check for numerically close roots.
    
    # Sort roots to make checking for duplicates easier
    sorted_roots = np.sort_complex(final_roots)
    is_repeated = np.zeros(n_final, dtype=bool)
    if n_final > 1:
        for i in range(1, n_final):
            if np.isclose(sorted_roots[i], sorted_roots[i-1]):
                # If a root is close to its predecessor, mark both as part of a repeated group
                is_repeated[i] = True
                is_repeated[i-1] = True

    max_norm_contrib = 0.0
    for i in range(n_final):
        root = sorted_roots[i]
        norm_contrib = np.abs(root)
        if is_repeated[i]:
            norm_contrib += 1.0  # Add 1 to norm for repeated roots
            
        if norm_contrib > max_norm_contrib:
            max_norm_contrib = norm_contrib

    infinity_norm = max_norm_contrib
    result = n_final * infinity_norm

    # Step 4: Output the final answer
    print(f"Identified the smallest n such that f(n) > 10.")
    print(f"This occurs at n = {n_final}, with f({n_final}) = {final_f_n:.4f}.")
    print(f"The infinity norm ||W_{n_final}||_infinity is calculated as {infinity_norm:.4f}.")
    print("\nFinal calculation:")
    print(f"{n_final} * ||W_{n_final}||_inf = {n_final} * {infinity_norm:.4f} = {result:.4f}")
    
    # Final answer in the requested format
    print(f"\n<<<{result:.4f}>>>")

solve_and_print()