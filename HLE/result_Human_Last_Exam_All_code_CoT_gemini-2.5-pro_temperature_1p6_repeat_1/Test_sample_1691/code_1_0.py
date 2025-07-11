import math

def get_integral_approximation():
    """
    This function computes the parameters for the analytical approximation of the given integral
    for the small epsilon regime.
    The approximation is of the form I(epsilon) ≈ C * epsilon^(-a).
    """
    
    # From the denominator h(x) = 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0,
    # we identify the leading term for small x, which is h(x) ≈ c0 * x^p.
    c0 = 9.0
    p = 5.0

    # The exponent 'a' in the approximation I(epsilon) ∝ epsilon^(-a)
    # is given by a = 1 - 1/p.
    a = 1.0 - (1.0 / p)

    # The coefficient C is given by the formula:
    # C = pi / (p * c0^(1/p) * sin(pi/p))
    C = math.pi / (p * math.pow(c0, 1.0/p) * math.sin(math.pi/p))

    # Output the results in a clear format.
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print("I(epsilon) ≈ C * epsilon^(-a)")
    print("\nBased on the leading term h(x) ≈ 9.0 * x^5.0:")
    
    # Output each number in the final equation as requested.
    # The numbers are the coefficient C and the exponent -a.
    print("\nThe parameters are:")
    print(f"Coefficient C = {C:.6f}")
    print(f"Exponent value -a = -{a}")
    
    print("\nThe final equation is:")
    print(f"I(epsilon) ≈ {C:.6f} * epsilon^(-{a})")


if __name__ == '__main__':
    get_integral_approximation()