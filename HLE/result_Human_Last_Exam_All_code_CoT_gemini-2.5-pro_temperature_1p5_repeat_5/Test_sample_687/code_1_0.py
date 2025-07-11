import spherogram
import sympy

def find_seifert_circle_lower_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    of the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # Step 1: Get the Link object for the 9_23 knot.
        knot = spherogram.Link('9_23')
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure the 'spherogram' and 'sympy' libraries are installed ('pip install spherogram sympy')")
        return

    # Step 2: Calculate the HOMFLY polynomial, P(l, m).
    # Spherogram uses the skein relation: l*P(+) + l^{-1}*P(-) + m*P(0) = 0.
    # The span is calculated with respect to the 'l' variable.
    p_homfly = knot.homfly_polynomial()

    # Identify the 'l' variable from the polynomial's symbols.
    variables = p_homfly.free_symbols
    l_var = None
    for v in variables:
        if v.name.lower() == 'l':
            l_var = v
    
    if l_var is None:
        print("Error: Could not identify the 'l' variable in the polynomial.")
        return

    print(f"The HOMFLY polynomial for the 9_23 knot is: P(l, m) = {p_homfly}\n")

    # Step 3: Calculate the span of the polynomial in the variable 'l'.
    # We convert the expression to a Poly object to easily get the degrees.
    p_poly_in_l = sympy.poly(p_homfly, l_var)
    
    # The terms method gives tuples of (powers, coefficient). We extract the powers.
    # For a multivariate poly, powers is a tuple, but since we specified one gen (l_var), it's a tuple of length 1.
    degrees = [term[0][0] for term in p_poly_in_l.terms()]
    
    max_degree = max(degrees)
    min_degree = min(degrees)
    
    print(f"The powers of the variable '{l_var.name}' in the polynomial are: {sorted(degrees)}")
    print(f"Maximum power of '{l_var.name}': {max_degree}")
    print(f"Minimum power of '{l_var.name}': {min_degree}\n")
    
    span_l = max_degree - min_degree
    print(f"The span of the polynomial is the difference between the max and min powers:")
    print(f"span_l(P) = {max_degree} - ({min_degree}) = {span_l}\n")
    
    # Step 4: Apply the inequality to find the lower bound for the number of Seifert circles, s.
    # The inequality is s >= span_l(P) / 2 + 1.
    lower_bound_float = span_l / 2 + 1
    lower_bound = int(lower_bound_float)

    print("The number of Seifert circles (s) is bounded by the inequality: s >= span_l(P) / 2 + 1")
    print("Substituting the numbers for the 9_23 knot:")
    print(f"s >= {span_l} / 2 + 1")
    print(f"s >= {span_l / 2} + 1")
    print(f"s >= {lower_bound}\n")
    
    print(f"Thus, a lower bound for the minimum number of Seifert circles is {lower_bound}.")

if __name__ == "__main__":
    find_seifert_circle_lower_bound()