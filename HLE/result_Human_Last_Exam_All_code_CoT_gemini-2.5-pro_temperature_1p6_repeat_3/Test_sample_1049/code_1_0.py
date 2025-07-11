import sympy

def find_closed_form():
    """
    Finds a closed form for the sum:
    Sum_{k=0 to n} (2k+1)**5 * C(2k, k) * C(2n-2k, n-k)
    using the generating functions method.
    """
    x, n = sympy.symbols('x n')

    # Define the generating function for C(2k, k)
    B0 = (1 - 4*x)**(-sympy.Rational(1, 2))

    # Define the operator (2*theta + 1), where theta = x*d/dx
    def apply_operator(f):
        return sympy.simplify(2 * x * sympy.diff(f, x) + f)

    # Iteratively apply the operator to find A_5(x)
    A = B0
    print("Deriving generating functions A_m(x):")
    print(f"A_0(x) = {A}")
    for m in range(5):
        A = apply_operator(A)
        print(f"A_{m+1}(x) = {A}")

    A5 = A

    # The generating function for the sum S_n is A_5(x) * B_0(x)
    S_series = sympy.simplify(A5 * B0)
    print("\nGenerating function for the sum S_n:")
    print(f"S(x) = A_5(x) * B_0(x) = {S_series}")

    # The resulting series is of the form P(x) / (1-4x)^6
    # Extract the numerator polynomial P(x)
    num, den = sympy.fraction(sympy.cancel(S_series))
    P5x = sympy.expand(num)
    
    print("\nNumerator polynomial P(x) in S(x) = P(x) / (1-4x)^6:")
    print(f"P(x) = {P5x}")

    # Get coefficients of P(x)
    p_coeffs = sympy.Poly(P5x, x).all_coeffs()
    p_coeffs.reverse()  # Start with the coefficient of x^0
    print(f"\nCoefficients of P(x) = c_0 + c_1*x + ...: {p_coeffs}")

    # The coefficient of x^k in (1-4x)^-m is C(k+m-1, m-1) * 4^k
    # Here, m=6. The coefficient of x^k in (1-4x)^-6 is C(k+5, 5) * 4^k
    # S_n = [x^n] S(x) = [x^n] (sum_{i=0}^{4} c_i * x^i) * (1-4x)^-6
    # S_n = sum_{i=0}^{4} c_i * [x^{n-i}] (1-4x)^-6
    # S_n = sum_{i=0}^{4} c_i * C(n-i+5, 5) * 4^(n-i)

    # Let's write the formula for S_n
    s_n_terms = []
    for i, c in enumerate(p_coeffs):
        term = c * sympy.binomial(n - i + 5, 5) * 4**(n - i)
        s_n_terms.append(term)
    
    Sn = sum(s_n_terms)
    
    # Simplify the expression
    Sn_simplified = sympy.simplify(Sn)
    
    # Extract the polynomial part of S_n = 4^n * Q(n)
    q_n, _ = sympy.apart(Sn_simplified / (4**n), n).as_numer_denom()
    q_n = sympy.expand(q_n)
    
    # The final closed form is 4^n * Q(n). Let's print the parts.
    
    print("\nThe closed form is S_n = 4**n * ( ... )")
    
    final_equation_poly_part = f""
    den = sympy.fraction(q_n)[1]
    num_poly = sympy.Poly(sympy.fraction(q_n)[0],n)
    
    # To match the requested format: "output each number in the final equation"
    # We will build the string for the polynomial Q(n).
    
    if den != 1:
        final_equation_poly_part += "("

    first_term = True
    for coeff, monom in zip(num_poly.coeffs(), num_poly.monoms()):
        if not first_term:
            if coeff > 0:
                final_equation_poly_part += " + "
            else:
                final_equation_poly_part += " - "
                coeff = -coeff
        else:
            if coeff < 0:
                final_equation_poly_part += "-"
                coeff = -coeff
        
        first_term = False
        
        final_equation_poly_part += f"{coeff}*n**{monom[0]}"

    if den != 1:
        final_equation_poly_part += f") / {den}"

    print(f"\nThe polynomial Q(n) is: {q_n}")
    print("\nFinal closed form expression:")
    print(f"S_n = 4**n * ({q_n})")

find_closed_form()