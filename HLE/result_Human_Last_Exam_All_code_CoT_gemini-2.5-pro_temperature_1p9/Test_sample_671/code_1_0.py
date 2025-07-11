import sympy

def get_fn_polynomial(n_val):
    """
    Generates the polynomial f_n(p) for a given integer n using the derived formula.

    Args:
        n_val (int): The index of the polynomial to generate. Must be >= 1.

    Returns:
        sympy.Expr: The sympy expression for the polynomial f_n(p).
    """
    if not isinstance(n_val, int) or n_val < 1:
        raise ValueError("n must be an integer greater than or equal to 1.")
        
    p = sympy.Symbol('p')
    
    # The closed-form formula for f_n(p) is (p^n - (1-p)^n) / (2*p - 1).
    # For the special case p=1/2, sympy's simplify/expand correctly computes the limit.
    formula = (p**n_val - (1 - p)**n_val) / (2 * p - 1)
    
    # Expand the formula to get the polynomial form.
    poly = sympy.expand(formula)
    
    return poly

def main():
    """
    Main function to demonstrate the formula by generating f_9(p)
    and printing its equation and coefficients.
    """
    n = 9
    
    # Generate the polynomial using the function
    fn_poly = get_fn_polynomial(n)
    
    # Create a Poly object to easily extract coefficients
    poly_obj = sympy.Poly(fn_poly, sympy.Symbol('p'))
    coefficients = poly_obj.all_coeffs()
    
    # Create the string representation of the final equation
    p_str = sympy.Symbol('p')
    equation_str = sympy.Eq(sympy.Function('f')(n, p_str), fn_poly)

    print(f"Based on the derived formula, the final equation for n={n} is:")
    # Using pretty print for a nicely formatted equation
    sympy.pprint(equation_str, use_unicode=True)

    print("\nAs requested, here is each number (coefficient) in the final equation, from the highest power of p to the constant term:")
    for coeff in coefficients:
        print(coeff)

if __name__ == "__main__":
    main()
