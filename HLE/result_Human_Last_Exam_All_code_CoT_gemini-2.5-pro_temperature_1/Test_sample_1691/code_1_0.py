import numpy as np

def solve():
    """
    This function calculates the coefficients for the asymptotic expansion of the integral I(epsilon)
    and prints the final analytical formula.
    """
    # Parameters from the integral denominator
    # g(x) = a*x^n1 + b*x^n2 + ...
    a = 9.0
    n1 = 5.0
    b = 5.0
    n2 = 6.0

    # --- Calculation of the first term: C1 * epsilon^(-(n1-1)/n1) ---
    # The first term comes from integrating 1 / (epsilon + a*x^n1)
    # The exponent of epsilon is -(n1-1)/n1
    p1_num = n1 - 1.0
    p1_den = n1
    p1 = p1_num / p1_den

    # The coefficient C1 is given by pi / (n1 * a^(1/n1) * sin(pi/n1))
    C1 = np.pi / (n1 * (a**(1.0/n1)) * np.sin(np.pi/n1))

    # --- Calculation of the second term: C2 * epsilon^(-(2*n1-n2-1)/n1) ---
    # The second term comes from integrating -b*x^n2 / (epsilon + a*x^n1)^2
    # The exponent of epsilon is -(2*n1 - n2 - 1)/n1
    p2_num = (2.0*n1 - n2 - 1.0)
    p2_den = n1
    p2 = p2_num / p2_den

    # The coefficient C2 is derived from the integral evaluation.
    # It evaluates to - (2*pi*b) / (n1^2 * a^((n2+1)/n1) * sin(pi*(n2-n1+1)/n1))
    # For our specific values, this simplifies to:
    C2 = (-2.0 * np.pi * b) / (n1**2 * a**((n2 + 1.0) / n1) * np.sin(np.pi * (n2 - n1 + 1.0) / n1))
    
    # Alternatively, using the simplified formula from the thought process:
    # C2 = -2.0 * np.pi / (5.0 * (9.0**(7.0/5.0)) * np.sin(2.0*np.pi/5.0))

    # --- Print the resulting formula ---
    print("An analytical formula that approximates I(epsilon) for the small epsilon regime is an asymptotic expansion.")
    print("The first two terms of this expansion are:")
    print(f"\nI(epsilon) \u2248 {C1:.4f} * \u03B5^(-{p1_num}/{p1_den}) + ({C2:.4f}) * \u03B5^(-{p2_num}/{p2_den})\n")
    print("This can also be written with decimal exponents as:")
    print(f"I(epsilon) \u2248 {C1:.4f} * \u03B5^(-{p1:.2f}) {C2:+.4f} * \u03B5^(-{p2:.2f})")

solve()