import spherogram
import sympy

def solve_knot_problem():
    """
    This function solves the knot theory problem by:
    1. Defining the braid from the problem description.
    2. Computing its closure to form a link.
    3. Identifying the components of the link.
    4. Calculating the Alexander polynomial for each component to determine its knot type.
    5. Printing the final result.
    """
    # The variable for the Alexander polynomial
    t = sympy.var('t')

    # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^{-1} in B_5.
    # In spherogram's notation, sigma_i is represented by the integer i,
    # and its inverse sigma_i^{-1} is represented by -i.
    braid_word = [1, 1, 2, 2, 3, -4]
    num_strands = 5

    print(f"Analyzing the closure of the braid on {num_strands} strands given by the word: {braid_word}")
    
    # Create the braid object
    braid = spherogram.Braid(num_strands, braid_word)

    # Compute the closure of the braid to get the link
    link = braid.closure()

    # Get the connected components of the link
    components = link.link_components
    print(f"The resulting link has {len(components)} components.")

    # Define known Alexander polynomials for comparison.
    # Note: These are standard forms. The computed ones might differ by a factor of +/-t^k.
    knot_map = {
        "Unknot": sympy.Integer(1),
        "Trefoil": t - 1 + t**-1,
        "Figure-8": -t + 3 - t**-1,
        "$5_1$ knot": t**2 - t + 1 - t**-1 + t**-2
    }

    def normalize_poly(p):
        """
        Normalizes a Laurent polynomial to a canonical form for comparison.
        The Alexander polynomial is defined up to multiplication by +/-t^k.
        This function converts it to a regular polynomial with a positive leading
        coefficient and lowest possible degree. It also handles the P(t) vs P(1/t) ambiguity.
        """
        if not isinstance(p, sympy.Expr) or p == 0:
            return sympy.Integer(0) if p == 0 else p

        # Ensure p is a poly in t
        if not p.has(t):
            return abs(p)

        p_poly = sympy.poly(p, t)
        # Make it a regular polynomial (no negative powers)
        lowest_degree = p_poly.lowest_deg()
        p_norm = sympy.expand(p_poly.as_expr() * t**(-lowest_degree))
        
        # Make leading coefficient positive
        if p_norm.as_poly(t).LC() < 0:
            p_norm = -p_norm
            
        # The Alexander polynomial P(t) is often symmetric, P(t) = P(1/t).
        # We choose the canonical form to be the one with the "smaller" representation.
        p_inv = sympy.expand(p_norm.subs(t, t**-1))
        p_inv_poly = sympy.poly(p_inv, t)
        p_inv_norm = sympy.expand(p_inv_poly.as_expr() * t**(-p_inv_poly.lowest_deg()))

        if str(p_inv_norm) < str(p_norm):
            return p_inv_norm
        else:
            return p_norm

    # Normalize the reference polynomials
    normalized_knot_map = {name: normalize_poly(poly) for name, poly in knot_map.items()}
    
    print("\nIdentifying each component by its Alexander polynomial:")
    
    unknot_count = 0
    other_component_name = None

    for i, component in enumerate(components):
        poly = component.alexander_polynomial()
        norm_poly = normalize_poly(poly)
        
        found_match = False
        for name, ref_poly in normalized_knot_map.items():
            if norm_poly == ref_poly:
                print(f"Component {i+1} is identified as: {name}")
                if name == "Unknot":
                    unknot_count += 1
                else:
                    other_component_name = name
                found_match = True
                break
        if not found_match:
            print(f"Component {i+1} is an unknown knot. Normalized polynomial: {norm_poly}")

    print("\n--- Conclusion ---")
    if unknot_count == 2 and other_component_name is not None:
        print(f"As stated in the problem, two components are unknots.")
        print(f"The other connected component is equivalent to the {other_component_name}.")
    elif len(components) == unknot_count:
         print("All components were found to be unknots.")
         other_component_name = "Unknot"
    else:
        print("The result does not match the problem's premise.")

solve_knot_problem()