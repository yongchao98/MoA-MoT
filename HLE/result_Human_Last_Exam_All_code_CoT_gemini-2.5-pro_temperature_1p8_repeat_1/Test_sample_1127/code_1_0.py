import sympy

def solve():
    """
    This function derives the minimal polynomial for the connective constant of the specified graph G.
    It uses the transfer matrix method for self-avoiding walks on a strip.
    """
    # z is the symbolic variable for the generating function, weighting each edge by its length.
    z = sympy.Symbol('z')

    # Define the transfer matrix T(z) based on enumerating SAW paths in a K4-e cell.
    # u0=(x,0), u1=(x,1) are inputs. v0=(x+1,0), v1=(x+1,1) are outputs.
    # T00: u0 -> v0 paths: z (len 1) + z^2 (len 2) + z^3 (len 3)
    T00 = z + z**2 + z**3
    # T01: u0 -> v1 paths: 2*z^2 (len 2) + 2*z^3 (len 3)
    T01 = 2*z**2 + 2*z**3
    # T10: u1 -> v0 paths: z (len 1) + 2*z^2 (len 2)
    T10 = z + 2*z**2
    # T11: u1 -> v1 paths: z (len 1) + z^2 (len 2) + z^3 (len 3)
    T11 = z + z**2 + z**3
    
    T = sympy.Matrix([[T00, T01], [T10, T11]])

    # The connective constant is found when the largest eigenvalue of T(z) is 1.
    # Find eigenvalues of T. λ is the eigenvalue.
    lam = sympy.Symbol('lam')
    char_poly = sympy.det(T - lam * sympy.eye(2))
    
    # Eigenvalues are the roots of the characteristic polynomial.
    # For a 2x2 matrix [[a,b],[c,d]], eigenvalues are (a+d)/2 ± sqrt(((a-d)/2)^2 + bc)
    # λ = (T00+T11)/2 ± sqrt(((T00-T11)/2)^2 + T01*T10)
    # Here T00 = T11, so λ = T00 ± sqrt(T01*T10)
    eigenvals = [T00 + sympy.sqrt(T01*T10), T00 - sympy.sqrt(T01*T10)]
    
    # We are interested in the largest eigenvalue for positive z, which corresponds to the '+' sign.
    lambda_max = eigenvals[0]
    
    # Solve lambda_max(z) = 1 for z. The solution z_c gives the connective constant mu = 1/z_c.
    # 1 = T00 + sqrt(T01*T10)
    # 1 - T00 = sqrt(T01*T10)
    # Square both sides to get a polynomial equation.
    lhs = (1 - T00)**2
    rhs = T01 * T10
    
    poly_z_eq = sympy.expand(lhs - rhs)

    # We need the polynomial for mu. Substitute z = 1/mu. Let mu=x.
    x = sympy.Symbol('x')
    poly_x_eq, _ = sympy.fraction(sympy.simplify(poly_z_eq.subs(z, 1/x)))
    
    # The result is the numerator. We make the leading coefficient positive.
    final_poly = sympy.poly(poly_x_eq, x)
    if final_poly.LC() < 0:
        final_poly = -final_poly

    # Print the equation term by term
    print("The minimal polynomial P(x) for the connective constant x=μ is P(x) = 0, where:")
    coeffs = final_poly.all_coeffs()
    degree = final_poly.degree()
    poly_str = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if coeff == 0:
            continue
        
        sign = "+ " if coeff > 0 else "- "
        coeff = abs(coeff)

        if i == 0:
            sign = "" if coeff > 0 else "- "

        term = ""
        if coeff == 1 and power != 0:
            term = f"x^{power}"
        else:
            term = f"{coeff}*x^{power}" if power > 1 else f"{coeff}*x" if power == 1 else str(coeff)
        
        poly_str.append(sign + term)

    print("P(x) = " + " ".join(poly_str))
    
    # Final output requirement: print each coefficient of the final polynomial.
    print("\nThe coefficients of the polynomial, from the highest degree term to the constant term, are:")
    # A degree 6 polynomial should have 7 coefficients. Some might be zero.
    all_c = final_poly.all_coeffs()
    # The problem asks to output the equation too
    print(f"{all_c[0]}*x^6 + ({all_c[1]})*x^5 + ({all_c[2]})*x^4 + ({all_c[3]})*x^3 + ({all_c[4]})*x^2 + ({all_c[5]})*x + ({all_c[6]}) = 0")


solve()