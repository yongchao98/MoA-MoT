import sympy

def solve_integral():
    """
    Constructs the symbolic analytical value of the integral.
    """
    # Define symbols
    s, a, b, i = sympy.symbols('s a b i', real=True)
    i = sympy.I  # Imaginary unit
    pi = sympy.pi
    zeta = sympy.zeta # Hurwitz Zeta function
    
    # Parameters from the derivation
    s_val = sympy.Rational(3, 2)
    a_val = sympy.Rational(3, 2)
    b_val = sympy.Rational(2, 3)

    # The complex arguments for the Hurwitz Zeta function
    q1 = a_val + i * b_val
    q2 = a_val - i * b_val
    
    # The analytical expression for the real part of the integral (J)
    J = sympy.sqrt(pi) * (zeta(s_val, q1) + zeta(s_val, q2))
    
    # The final integral I = iJ
    I = i * J

    # Format the equation string
    # We want to show the numbers that make up the final expression.
    J_str_expr = f"sqrt(pi) * (zeta({s_val}, {a_val} + {b_val}*I) + zeta({s_val}, {a_val} - {b_val}*I))"
    
    print(f"The analytical value of the integral is I = i * J where:")
    print(f"J = {J_str_expr}")
    # We can also print the full expression for I
    I_str_expr = f"I = I * {J_str_expr}"
    print("\nSymbolically, the integral is:")
    print(I_str_expr)


if __name__ == "__main__":
    solve_integral()
