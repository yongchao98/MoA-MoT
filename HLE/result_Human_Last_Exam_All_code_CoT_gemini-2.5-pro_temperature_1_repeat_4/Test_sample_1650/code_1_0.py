import sympy

def solve_overlap_integral():
    """
    Calculates the overlap integral for two 2s orbitals in H2+
    using symbolic mathematics and prints the final analytical expression.
    """
    # Define symbols
    # w is a dimensionless variable representing zeta * R
    w = sympy.Symbol('w', positive=True)
    # lambda is the elliptical coordinate
    lam = sympy.Symbol('lambda')

    # The integration over phi and mu results in a polynomial in lambda, I_mu(lam).
    # The derivation shows this polynomial is:
    # (w^2/2)λ^4 - 4wλ^3 + (8 - w^2/3)λ^2 + (4w/3)λ + (w^2/10 - 8/3)
    poly_coeffs = [
        (w**2/10 - sympy.S(8)/3),  # Coefficient for lambda^0
        (sympy.S(4)/3 * w),           # Coefficient for lambda^1
        (8 - w**2/3),                 # Coefficient for lambda^2
        -4*w,                         # Coefficient for lambda^3
        w**2/2                        # Coefficient for lambda^4
    ]

    # The next step is to integrate each term of the polynomial multiplied
    # by exp(-w*lambda/2) from lambda = 1 to infinity.
    # We define a function for the auxiliary integral A_n(a) = integral(x^n * exp(-a*x))
    a = w / 2
    def auxiliary_integral(n, var_a):
        """Computes the definite integral of lambda^n * exp(-a*lambda) from 1 to oo."""
        return sympy.integrate(lam**n * sympy.exp(-var_a * lam), (lam, 1, sympy.oo))

    # Calculate the total value of the lambda integral
    total_lambda_integral = sum(
        c * auxiliary_integral(n, a) for n, c in enumerate(poly_coeffs)
    )

    # The overlap integral S is obtained by multiplying by the remaining constants
    # from the normalization and the coordinate transformation.
    # This factor is (w^3 / 128).
    S_expression = (w**3 / 128) * total_lambda_integral

    # Simplify the final expression
    S_final = sympy.expand(sympy.simplify(S_expression))

    # Extract the polynomial part P(w) where S = exp(-w/2) * P(w)
    poly_part = sympy.simplify(S_final / sympy.exp(-w/2))

    # Format the output to clearly show the final equation
    print("The overlap integral S for two 2s orbitals is given by:")
    print("S(w) = exp(-w/2) * P(w)")
    print("where w = ζR (ζ is effective nuclear charge, R is internuclear distance)")
    print("\nThe polynomial P(w) is:")
    
    # We want to print the polynomial in a clean, readable format.
    # Example: 1 + 1/2*w + 1/12*w**2 - 1/120*w**3 - 1/120*w**4
    final_poly = sympy.Poly(poly_part, w)
    terms = []
    for i in range(final_poly.degree(), -1, -1):
        coeff = final_poly.coeff_monomial(w**i)
        if coeff == 0:
            continue
        
        sign = "+" if coeff > 0 else "-"
        # Don't print + for the first term
        if not terms and sign == "+":
            sign = ""
            
        coeff_abs = abs(coeff)
        
        if i == 0:
            term_str = f"{coeff_abs}"
        elif i == 1:
            if coeff_abs == 1:
                 term_str = "w"
            else:
                 term_str = f"{coeff_abs}*w"
        else:
            if coeff_abs == 1:
                 term_str = f"w**{i}"
            else:
                 term_str = f"{coeff_abs}*w**{i}"

        terms.append(f"{sign} {term_str}")

    print("P(w) =", " ".join(terms))
    print("\nSo the full equation for the overlap integral is:")
    print(f"S(w) = exp(-w/2) * ({poly_part})")


solve_overlap_integral()