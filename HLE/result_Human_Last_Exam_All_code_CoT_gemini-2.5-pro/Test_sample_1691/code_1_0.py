import math

def solve_integral_approximation():
    """
    This function calculates the analytical approximation for the integral
    I(eps) = integral_0^15 dx / (eps + 9*x^5 + 5*x^6 + 9*x^8)
    for the small epsilon regime.
    """

    # From the problem, the dominant term in the polynomial for small x is 9*x^5.
    # So, we have f(x) ~= c * x^p
    p = 5.0
    c = 9.0

    # The approximation for the integral is of the form I(eps) ~= C * eps^k.

    # 1. Calculate the exponent k
    # The derivation shows that k = 1/p - 1 = (1 - p) / p.
    exponent = (1.0 - p) / p

    # 2. Calculate the coefficient C
    # The derivation shows C = c^(-1/p) * integral_0^inf(1/(1+y^p) dy).
    # The integral part is equal to (pi/p) / sin(pi/p).
    integral_constant = (math.pi / p) / math.sin(math.pi / p)
    coefficient = (c**(-1.0 / p)) * integral_constant

    # 3. Print the final analytical formula with the calculated numbers.
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(epsilon) ~= {coefficient} * epsilon^({exponent})")
    print("\nWhere the numbers in the equation are:")
    print(f"Coefficient: {coefficient}")
    print(f"Exponent: {exponent}")

solve_integral_approximation()