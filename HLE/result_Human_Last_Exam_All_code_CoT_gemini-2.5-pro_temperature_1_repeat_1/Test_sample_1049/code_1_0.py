import sympy

def main():
    """
    This script computes the closed form for the given sum using sympy
    and prints the result as a formatted equation.
    """
    n = sympy.Symbol('n')

    # The expression for Q(n) is derived from the generating function method,
    # where S_n = 4^n * Q(n).
    Q_n_expr = (16 * sympy.binomial(n + 1, 5) +
                296 * sympy.binomial(n + 2, 5) +
                516 * sympy.binomial(n + 3, 5) +
                116 * sympy.binomial(n + 4, 5) +
                sympy.binomial(n + 5, 5))

    # Expand the expression to get the polynomial form
    poly_Q_n = sympy.expand(Q_n_expr)
    
    # The polynomial Q(n) has a factor of (n+1).
    # We find the other polynomial factor, P(n) = Q(n)/(n+1).
    p_n = sympy.expand(poly_Q_n / (n + 1))
    
    # To present the result with integer coefficients, we find the least common
    # multiple of the denominators of the coefficients of P(n).
    p = sympy.Poly(p_n, n)
    coeffs = p.all_coeffs()
    common_denom = sympy.lcm([sympy.Rational(c).q for c in coeffs])
    
    # Scale P(n) by the common denominator to get a polynomial with integer coefficients.
    scaled_p_n = sympy.expand(p_n * common_denom)
    
    p_scaled = sympy.Poly(scaled_p_n, n)
    coeffs_scaled = p_scaled.all_coeffs()
    degree_scaled = p_scaled.degree()
    
    # Build the string for the polynomial with integer coefficients.
    inner_poly_str_scaled = ""
    for i, c in enumerate(coeffs_scaled):
        power = degree_scaled - i
        if c == 0:
            continue
        
        # Add the sign (+ or -)
        if i > 0:
            if c > 0:
                inner_poly_str_scaled += " + "
            else:
                inner_poly_str_scaled += " - "
                c = abs(c)
        elif c < 0:
            inner_poly_str_scaled += "-"
            c = abs(c)

        # Add the integer coefficient
        inner_poly_str_scaled += str(int(c))
        
        # Add the variable n with its power
        if power > 0:
            inner_poly_str_scaled += "*n"
            if power > 1:
                inner_poly_str_scaled += f"^{power}"

    # Construct the final equation string.
    equation = f"S_n = 4^n * (n + 1) * ({inner_poly_str_scaled}) / {common_denom}"
    print(equation)

main()