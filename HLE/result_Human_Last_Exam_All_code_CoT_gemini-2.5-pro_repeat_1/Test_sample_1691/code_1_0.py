import math

def solve_integral_approximation():
    """
    This function calculates the constants for the asymptotic formula of the integral
    I(epsilon) = integral_0^15 (1 / (epsilon + 9*x^5 + 5*x^6 + 9*x^8)) dx
    The formula is of the form I(epsilon) ~ C1 * epsilon^C2 for small epsilon.
    """
    # The dominant term in the denominator for x -> 0 is g(x) ~ a * x^p
    a = 9.0
    p = 5.0

    # The exponent C2 is given by the scaling law
    C2 = (1.0 - p) / p

    # The coefficient C1 is derived from the substitution and a standard integral
    # C1 = a^(-1/p) * integral_0^inf(1/(1+u^p)) du
    # The standard integral is equal to (pi/p) / sin(pi/p)
    integral_part = (math.pi / p) / math.sin(math.pi / p)
    C1 = math.pow(a, -1.0 / p) * integral_part

    # Print the derived analytical formula with the computed constants.
    print(f"The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(epsilon) = C1 * epsilon^C2")
    print(f"where the numbers are:")
    print(f"C1 = {C1}")
    print(f"C2 = {C2}")
    print(f"Thus, the final equation is:")
    print(f"I(epsilon) = {C1} * epsilon^({C2})")

solve_integral_approximation()