import math

def solve():
    """
    This function calculates the coefficients for the analytical approximation of the integral.
    """

    # From the problem, the integral is I(epsilon) = integral_0^15 of 1 / (epsilon + 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0) dx.
    # For small x, the polynomial is dominated by the term with the smallest power of x.
    # The approximated denominator is epsilon + a * x^n.
    a = 9.0
    n = 5.0

    # The integral for small epsilon can be approximated by a formula of the form:
    # I(epsilon) ~= C * epsilon^p

    # The exponent p is given by the formula (1-n)/n.
    p = (1.0 - n) / n

    # The coefficient C is derived from the analytical solution of the approximated integral:
    # C = (pi / (n * sin(pi/n))) * a^(-1/n)
    C = (math.pi / (n * math.sin(math.pi / n))) * (a**(-1.0 / n))

    # The final analytical formula is I(epsilon) ~= C * epsilon^p.
    # We print the formula with the calculated numerical values for C and p.
    print(f"The approximated analytical formula is:")
    print(f"I(epsilon) ~= {C:.4f} * epsilon^({p:.1f})")
    print("\nwhere the numbers in the equation are:")
    print(f"Coefficient C = {C:.8f}")
    print(f"Exponent p = {p:.1f}")

solve()