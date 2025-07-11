import math

def solve_integral_approximation():
    """
    This function calculates the coefficients for the analytical approximation of the integral
    I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9x^5 + 5x^6 + 9x^8) dx
    for small epsilon.
    The approximation is of the form I(epsilon) ≈ C * epsilon^b.
    """

    # The polynomial in the denominator is P(x) = 9.0*x^5 + 5.0*x^6 + 9.0*x^8.
    # For small x, the behavior is dominated by the term with the lowest power.
    # The coefficient of the lowest power term (x^p).
    c_p = 9.0
    # The lowest power p.
    p = 5.0

    # The exponent 'b' in the approximation is given by b = -(p - 1) / p.
    b = -(p - 1) / p

    # The coefficient 'C' is given by C = (1 / c_p^(1/p)) * ( (pi/p) / sin(pi/p) ).
    # This comes from the analytical solution of the approximated integral.
    constant_integral_part = (math.pi / p) / math.sin(math.pi / p)
    C = (1 / (c_p**(1/p))) * constant_integral_part

    # Output the results in a clear format.
    print("The analytical formula that approximates I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) ≈ C * epsilon^b")
    print("\nBased on the provided integral, the calculated values for the formula are:")
    print(f"C = {C}")
    print(f"b = {b}")
    print("\nSo, the final approximate formula is:")
    # We output each number that makes up the final equation as requested.
    print(f"I(epsilon) ≈ {C} * epsilon^({b})")

solve_integral_approximation()
