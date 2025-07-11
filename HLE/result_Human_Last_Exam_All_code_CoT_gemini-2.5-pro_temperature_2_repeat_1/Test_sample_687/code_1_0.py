import sympy
# We need the pyknotid library. If not installed, run: pip install pyknotid
try:
    import pyknotid.catalogue as cat
except ImportError:
    print("The 'pyknotid' library is required. Please install it by running: pip install pyknotid")
    exit()

def solve_seifert_circle_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles of the 9_23 knot
    using its HOMFLY polynomial.
    """
    # Step 1: Get the specified knot from the catalogue.
    try:
        knot = cat.Knot(9, 23)
    except Exception as e:
        print(f"Failed to load knot 9_23. Error: {e}")
        return

    # Step 2: Compute the HOMFLY polynomial for the knot.
    # The pyknotid library returns a sympy expression in variables 'a' and 'z'.
    homfly_poly = knot.homfly_polynomial()
    print(f"The HOMFLY polynomial P(a, z) for the {knot.name} knot is:")
    print(f"P(a, z) = {homfly_poly}\n")

    # Step 3: Find the span of the polynomial with respect to the variable 'a'.
    # This requires finding the maximum and minimum degree of 'a'.
    variables = homfly_poly.free_symbols
    try:
        a_var = next(v for v in variables if str(v) == 'a')
    except StopIteration:
        print("Could not find variable 'a' in the polynomial.")
        return
        
    p_a = sympy.Poly(homfly_poly, a_var)
    degrees_of_a = [monomial[0] for monomial in p_a.monoms()]
    
    max_degree = max(degrees_of_a)
    min_degree = min(degrees_of_a)
    span_a = max_degree - min_degree
    
    print("To find the lower bound, we inspect the powers of variable 'a':")
    print(f" - The maximum power of 'a' is: {max_degree}")
    print(f" - The minimum power of 'a' is: {min_degree}\n")
    
    print("The span of the polynomial in 'a' is the difference between these powers:")
    print(f"span_a(P) = {max_degree} - {min_degree} = {span_a}\n")

    # Step 4: Apply the Morton-Franks-Williams inequality to find the lower bound.
    # The inequality is s(K) >= span_a(P) / 2 + 1.
    bound_part_1 = span_a // 2
    final_bound = bound_part_1 + 1
    
    print("Using the Morton-Franks-Williams inequality, s(K) >= span_a(P)/2 + 1, we calculate the bound:")
    print(f"s(K) >= {span_a} / 2 + 1 = {bound_part_1} + 1 = {final_bound}\n")
    
    print(f"Thus, a lower bound for the minimum number of Seifert circles for the {knot.name} knot is {final_bound}.")

if __name__ == '__main__':
    solve_seifert_circle_bound()