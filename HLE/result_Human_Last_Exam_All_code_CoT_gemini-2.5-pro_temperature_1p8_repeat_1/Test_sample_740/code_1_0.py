import sympy

def find_lambda_equation(n):
    """
    Finds the polynomial equation for lambda for which the integral equation has no solution,
    for a given integer n.

    The problem reduces to a linear system (I - lambda * M)c = f, where 1/lambda
    are the eigenvalues of the matrix M.
    Let mu be an eigenvalue of M. The characteristic equation is det(M - mu*I) = 0.
    Since mu = 1/lambda, we can find the equation for lambda.
    """
    if not isinstance(n, int) or n <= 0:
        print("n must be a positive integer.")
        return

    # In our derivation, we arrive at a matrix M where M_mk = 1/(n+k-m)
    # for m, k from 0 to n-1.
    M_entries = {}
    for m in range(n):
        for k in range(n):
            M_entries[(m, k)] = sympy.Rational(1, n + k - m)
            
    M = sympy.Matrix(n, n, lambda m, k: M_entries[(m, k)])

    # Let mu be the symbolic variable for the eigenvalues of M
    mu = sympy.Symbol('mu')
    
    # Calculate the characteristic polynomial det(M - mu*I)
    char_poly = M.charpoly(mu)
    
    # Let lambda_sym be the symbolic variable for lambda
    lambda_sym = sympy.Symbol('lambda')

    # Substitute mu = 1/lambda_sym into the characteristic polynomial
    lambda_eq_rational = char_poly.as_expr().subs(mu, 1/lambda_sym)
    
    # Convert the rational expression in lambda to a polynomial equation
    # We multiply by lambda_sym^n to clear the denominator
    lambda_eq_poly, _ = sympy.fraction(sympy.simplify(sympy.expand(lambda_eq_rational * lambda_sym**n)))
    
    # Ensure the leading coefficient is positive
    if lambda_eq_poly.LC() < 0:
        lambda_eq_poly = -lambda_eq_poly
        
    # Get the coefficients of the polynomial in lambda
    coeffs = lambda_eq_poly.all_coeffs()
    
    # Build the equation string
    equation_str = ""
    for i, coeff in enumerate(coeffs):
        power = n - i
        # Format the coefficient and variable part
        if power > 1:
            term_str = f"{coeff}*lambda^{power}"
        elif power == 1:
            term_str = f"{coeff}*lambda"
        else: # power == 0
            term_str = f"{coeff}"
        
        # Add sign
        if i == 0:
            equation_str += term_str
        else:
            if coeff > 0:
                equation_str += f" + {term_str}"
            else:
                # sympy coeff will already have a sign
                equation_str += f" - {-coeff}"
                if power > 1:
                    equation_str += f"*lambda^{power}"
                elif power == 1:
                    equation_str += f"*lambda"


    print(f"For n = {n}, the equation has no solution if lambda is a root of the following polynomial equation:")
    print(f"{equation_str} = 0")
    print("\nFor the specific final answer required, the equation for n=3 is:")
    # For the required output format, extracting the equation again for the final tag
    final_eq = sympy.Poly(lambda_eq_poly, lambda_sym)
    final_eq_str = str(final_eq.as_expr()) + " = 0"
    print(f"<<<{final_eq_str}>>>")


if __name__ == '__main__':
    # The problem is stated for a general n. We demonstrate for n=3.
    # To try other values, change the value of n below.
    n_value = 3
    find_lambda_equation(n_value)
