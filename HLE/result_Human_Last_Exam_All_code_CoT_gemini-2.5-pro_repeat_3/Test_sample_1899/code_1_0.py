import numpy as np

def solve():
    """
    Calculates the infimum and supremum for the given problem.
    """
    # Define the coefficients of the polynomial f(x) = x^3 + x^2 + x - 1
    f_coeffs = [1, 1, 1, -1]

    # Find the roots of the polynomial
    roots = np.roots(f_coeffs)

    # Identify the real root xi
    xi = None
    for r in roots:
        if np.isreal(r):
            xi = r.real
            break

    if xi is None:
        print("Real root not found.")
        return

    # f'(x) = 3x^2 + 2x + 1
    # The limit L is f'(xi) / 22
    f_prime_xi = 3 * xi**2 + 2 * xi + 1
    limit_val = f_prime_xi / 22

    # For this problem, the infimum and supremum are both equal to the limit.
    inf_val = limit_val
    sup_val = limit_val

    print(f"The root xi is approximately: {xi}")
    print(f"The value of f'(xi) is approximately: {f_prime_xi}")
    print(f"The limit L = f'(xi)/22 is approximately: {limit_val}")
    print(f"The infimum is: {inf_val}")
    print(f"The supremum is: {sup_val}")

solve()