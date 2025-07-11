import sympy

def solve_cohomology_dimension():
    """
    Calculates the dimension of the 9th cohomology group H^9(M, Q)
    based on the characteristic polynomial of the E6 quaternionic reflection arrangement.
    """
    # Define the symbolic variable for the polynomial
    t = sympy.Symbol('t')
    q = sympy.Symbol('q')

    # The characteristic polynomial of the E6 quaternionic arrangement in H^4
    # chi(q) = q^4 - 36q^3 + 450q^2 - 2040q + 3159
    chi_poly = q**4 - 36*q**3 + 450*q**2 - 2040*q + 3159
    
    # The argument for the characteristic polynomial in the formula for the Poincare series
    q_arg = (1 - t**3) / (1 - t**2)
    
    # Substitute the argument into the characteristic polynomial
    chi_substituted = chi_poly.subs(q, q_arg)
    
    # The formula for the Poincare polynomial of the affine complement M is:
    # P_M(t) = (1+t^3)/(1-t^2) * chi( (1-t^3)/(1-t^2) )
    poincare_poly_expr = ((1 + t**3) / (1 - t**2)) * chi_substituted
    
    # We need to find the coefficient of t^9 in the series expansion of P_M(t).
    # We expand the series up to a degree higher than 9, for instance, 10.
    series_expansion = sympy.series(poincare_poly_expr, t, 0, 10)
    
    # Extract the coefficient of t^9
    betti_9 = series_expansion.coeff(t, 9)
    
    # The problem asks to output the final equation.
    # The Poincare polynomial is the sum of Betti numbers times powers of t.
    # P(t) = b0 + b1*t + b2*t^2 + ...
    # We are looking for b9.
    
    # For clarity in the output, let's print the intermediate and final values.
    # The characteristic polynomial is:
    # q^4 - 36*q^3 + 450*q^2 - 2040*q + 3159
    # The formula for the Poincare polynomial is (1+t^3)/(1-t^2) * chi((1-t^3)/(1-t^2))
    # We are calculating the coefficient of t^9 in its series expansion.
    
    print("The characteristic polynomial is chi(q) = q^4 - 36*q^3 + 450*q^2 - 2040*q + 3159.")
    print("The Poincare polynomial P(t) is given by the formula (1+t^3)/(1-t^2) * chi((1-t^3)/(1-t^2)).")
    print("We need to find the coefficient of t^9 in the Taylor series of P(t).")
    print(f"The dimension of the ninth cohomology group H^9(M, Q) is the coefficient of t^9.")
    print(f"dim H^9(M, Q) = {betti_9}")

solve_cohomology_dimension()