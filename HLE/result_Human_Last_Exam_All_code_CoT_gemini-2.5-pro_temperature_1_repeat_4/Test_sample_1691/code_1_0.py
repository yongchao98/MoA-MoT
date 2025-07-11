import math

def solve_integral_approximation():
    """
    Calculates the analytical approximation for the integral I(epsilon) for small epsilon.
    The integral is I(eps) = integral from 0 to 15.00 of 1 / (eps + 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0) dx.
    For small epsilon, the approximation is of the form C * epsilon^k.
    This function calculates C and k and prints the resulting formula.
    """
    
    # Parameters from the integral's denominator:
    # f(x) = 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0
    
    # The dominant term for small x is the one with the lowest power of x.
    a = 9.0
    p = 5.0
    
    # The other terms in the polynomial.
    c1 = 5.0
    p1 = 6.0
    c2 = 9.0
    p2 = 8.0
    
    # The upper limit of integration.
    L = 15.00

    print("The integral to approximate is:")
    print(f"I(eps) = \u222B[0, {L}] dx / (eps + {a}*x^{p} + {c1}*x^{p1} + {c2}*x^{p2})")
    print("\nFor small epsilon, the leading-order approximation is obtained by considering the dominant behavior of the integrand near x=0.")
    print(f"This is determined by the lowest power of x in the denominator, which is the term: {a} * x^{p}")
    
    # The analytical approximation is of the form C * epsilon^k
    # Calculate the exponent k
    k = (1 / p) - 1

    # Calculate the coefficient C
    # The formula is C = (pi * a^(-1/p)) / (p * sin(pi/p))
    C = (math.pi * (a**(-1/p))) / (p * math.sin(math.pi / p))

    print("\nThe resulting analytical formula has the form: I(eps) \u2248 C * eps^k")
    print("\nPlugging in the values from the dominant term:")
    print(f"a = {a}")
    print(f"p = {p}")
    
    print("\nWe can calculate the exponent k and the coefficient C:")
    print(f"k = (1/{p}) - 1 = {k}")
    print(f"C = (\u03C0 * {a}^(-1/{p})) / ({p} * sin(\u03C0/{p}))")
    print(f"C \u2248 {C:.5f}")

    print("\nTherefore, the final analytical formula that approximates I(epsilon) is:")
    print(f"I(eps) \u2248 {C:.5f} * eps^({k})")

solve_integral_approximation()