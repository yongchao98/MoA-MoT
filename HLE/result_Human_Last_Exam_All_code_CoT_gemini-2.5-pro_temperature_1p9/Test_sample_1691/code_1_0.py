import numpy as np

def solve():
    """
    This function develops and prints an analytical approximation for the integral
    I(epsilon) for the small epsilon regime.
    """

    # The integral is dominated by the behavior near x=0.
    # The denominator is approximated by epsilon + c * x^p, where from the problem:
    c = 9.0
    p = 5.0

    # The approximate formula is I(epsilon) ~ K * epsilon^n.
    # The exponent n is given by 1/p - 1.
    n = 1.0/p - 1.0

    # The coefficient K is given by c**(-1/p) * integral_0^inf(1/(1+u**p))du
    # The integral evaluates to (pi/p)/sin(pi/p).
    K = (c**(-1.0/p)) * (np.pi/p) / np.sin(np.pi/p)

    print("The analytical formula that approximates I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) \u2248 K * \u03B5^n\n")

    print("The constants K and n are derived from the leading term of the polynomial in the denominator (9.0 * x^5.0).")
    
    # We output each number in the final equation by first stating what they represent.
    print(f"The coefficient from the dominant term is c = {c}")
    print(f"The power from the dominant term is p = {p}\n")
    
    print("This leads to the following values for K and n:")
    print(f"Coefficient K = (1 / ({c}^(1/{p}) * {p} * sin(\u03C0/{p}))) * \u03C0 \u2248 {K:.10f}")
    print(f"Exponent n = 1/{p} - 1 = {n}\n")

    print("The final approximate equation is:")
    print(f"I(\u03B5) \u2248 {K:.10f} * \u03B5^({n})")


solve()