import numpy as np

def solve_integral_approximation():
    """
    Calculates the parameters for the analytical approximation of the integral
    I(epsilon) = integral_0^15 (1 / (epsilon + 9x^5 + 5x^6 + 9x^8)) dx
    for the small epsilon regime.
    """

    # The approximation is based on the lowest power term in the polynomial, a*x^n.
    # From the problem, this is 9.0 * x^5.0.
    a = 9.0
    n = 5.0

    # The exponent 'p' in the approximation I(epsilon) ~ C * epsilon^p is given by
    # p = -(n - 1) / n
    p = -(n - 1.0) / n

    # The coefficient 'C' is given by C = a^(-1/n) * integral_0^inf(1 / (1 + u^n)) du.
    # The integral has a known solution: (pi/n) / sin(pi/n).
    integral_val = (np.pi / n) / np.sin(np.pi / n)
    C = np.power(a, -1.0 / n) * integral_val

    # Print the final analytical formula with the calculated coefficients.
    # The formula is I(epsilon) approx C * epsilon^p.
    # We output all the numbers that form the final equation.
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(\u03B5) \u2248 {C:.6f} * \u03B5^({p})")
    print("\nWhere the numbers in the equation are:")
    print(f"Coefficient C = {C}")
    print(f"Exponent p = {p}")

solve_integral_approximation()