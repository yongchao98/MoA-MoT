import sympy

def solve_for_lambda_equation(n):
    """
    This function finds the polynomial equation for lambda for a given n.
    The roots of this equation are the values of lambda for which the integral
    equation has no solution.

    Args:
        n (int): The exponent in the integral equation kernel.
    """
    if not isinstance(n, int) or n < 1:
        print("n must be a positive integer.")
        return

    # Define the matrix A
    A = sympy.zeros(n, n)
    for j in range(n):
        for k in range(n):
            A[j, k] = sympy.Rational(1, j + n - k)

    # Define the symbol for the eigenvalues of A
    mu = sympy.Symbol('mu')
    
    # Calculate the characteristic polynomial for mu: det(A - mu*I)
    I = sympy.eye(n)
    char_poly_mu = (A - mu * I).det()
    
    # The polynomial is often computed as det(mu*I - A).
    # The roots are the same. Let's use the standard form from sympy.
    char_poly_mu_sympy = A.charpoly(mu)

    # Get the coefficients of the polynomial in mu
    # p.all_coeffs() gives coeffs from highest power to lowest
    coeffs_mu = char_poly_mu_sympy.all_coeffs()
    
    # The equation for lambda is obtained by substituting mu = 1/lambda
    # and multiplying by lambda^n. This reverses the order of coefficients.
    coeffs_lambda = list(reversed(coeffs_mu))
    
    # Build and print the equation for lambda
    lambda_sym = sympy.Symbol('λ')
    equation = 0
    for i, coeff in enumerate(coeffs_lambda):
        equation += coeff * (lambda_sym ** i)
        
    print(f"For n = {n}, the matrix A is:")
    sympy.pprint(A)
    print("\nThe equation for λ is:")
    equation_str = sympy.pretty(sympy.Eq(equation, 0), use_unicode=True)
    print(equation_str)
    
    print("\nThe coefficients of the polynomial in λ (from highest power to lowest) are:")
    # The coefficients are already reversed, so we just need to reverse again for printing
    # from highest power to lowest. Or just print coeffs_lambda[::-1].
    final_coeffs = list(reversed(coeffs_lambda))
    for coeff in final_coeffs:
        print(coeff)
    
    print("\nTo express the equation more clearly, we can multiply by a common denominator.")
    # Use sympy.fraction to get numerator and denominator
    numer, denom = sympy.fraction(sympy.cancel(equation))
    final_equation = sympy.Eq(numer, 0)
    print("Simplified equation:")
    final_equation_str = sympy.pretty(final_equation, use_unicode=True)
    print(final_equation_str)

    print("\nThe numbers in the final simplified equation are the coefficients:")
    poly_final = sympy.Poly(final_equation.lhs, lambda_sym)
    for coeff in poly_final.all_coeffs():
        print(coeff)


# We will run the code for n=3 as a demonstration.
solve_for_lambda_equation(3)