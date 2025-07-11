import math

def solve_integral_approximation():
    """
    Calculates the parameters for the analytical approximation of the integral
    I(epsilon) for small epsilon.
    The integral is I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9*x^5 + 5*x^6 + 9*x^8) dx.
    The approximation has the form I(epsilon) ~ A * epsilon^B.
    """

    # Parameters from the denominator's leading term for small x: c1*x^p1
    c1 = 9.0
    p1 = 5.0

    # The exponent B is given by (1/p1) - 1
    exponent_B = (1.0 / p1) - 1.0

    # The coefficient A is given by (pi / sin(pi/p1)) / (p1 * c1^(1/p1))
    pi_over_sin_term = math.pi / math.sin(math.pi / p1)
    denominator_term = p1 * math.pow(c1, 1.0 / p1)
    coefficient_A = pi_over_sin_term / denominator_term
    
    # Print the resulting analytical formula
    print("The analytical formula that approximates I(epsilon) for the small epsilon regime is:")
    # The final equation string is constructed below
    # The numbers in the equation are coefficient_A and exponent_B
    final_equation = f"I(epsilon) ≈ {coefficient_A} * epsilon^({exponent_B})"
    print(final_equation)

if __name__ == "__main__":
    solve_integral_approximation()
