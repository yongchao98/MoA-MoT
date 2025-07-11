import sympy

def find_seifert_circle_lower_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles for the 9_23 knot
    using its HOMFLY polynomial.
    """
    
    # The HOMFLY polynomial for the 9_23 knot, P(v, z), from the KnotInfo database.
    # The inequality we use applies to this standard convention.
    homfly_poly_str = "v**-8 - v**-6 - v**-4*z**2 + 2*v**-4 + v**-2*z**2 - v**-2"
    
    print(f"The HOMFLY polynomial for the 9_23 knot is P(v, z) = {homfly_poly_str.replace('**','^')}")
    
    # Define symbolic variables and parse the polynomial string
    v, z = sympy.symbols('v z')
    poly = sympy.sympify(homfly_poly_str)
    
    # To find the span in v, we treat it as a polynomial in v
    poly_in_v = sympy.poly(poly, v)
    
    # Get all the powers of v that appear in the polynomial
    degrees_v = [monom[0] for monom in poly_in_v.monoms()]
    
    # Find the minimum and maximum degree of v
    min_degree = min(degrees_v)
    max_degree = max(degrees_v)
    
    # The span is the difference between the maximum and minimum degree
    span_v = max_degree - min_degree
    
    print(f"The powers of the variable 'v' in the polynomial are: {sorted(list(set(degrees_v)))}")
    print(f"The maximum power of v is: {max_degree}")
    print(f"The minimum power of v is: {min_degree}")
    print(f"The span of the polynomial in v is max_degree - min_degree = {max_degree} - ({min_degree}) = {span_v}")
    
    # The Cromwell-Lickorish-Millet inequality gives a lower bound for the number of Seifert circles s(K):
    # s(K) >= span_v(P)/2 + 1
    lower_bound = span_v / 2 + 1
    
    print("\nA lower bound for the minimum number of Seifert circles is given by the formula: span/2 + 1")
    print("Plugging in the numbers for the final equation:")
    print(f"{span_v} / 2 + 1 = {int(span_v / 2)} + 1 = {int(lower_bound)}")
    print(f"\nTherefore, a lower bound for the minimum number of Seifert circles of the 9_23 knot is {int(lower_bound)}.")

if __name__ == '__main__':
    find_seifert_circle_lower_bound()