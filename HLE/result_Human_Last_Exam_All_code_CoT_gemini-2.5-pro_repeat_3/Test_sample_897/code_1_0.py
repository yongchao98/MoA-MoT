import spherogram
import sympy

def solve_knot_problem():
    """
    This function solves the problem by calculating the two required knot theory quantities
    and finding their difference.
    
    Note: This script requires the 'spherogram' and 'sympy' libraries.
    You can install them using pip:
    pip install spherogram sympy
    """
    
    # --- Part 1: Analyze K1 = 10_74 ---
    
    # Get the knot K1 = 10_74
    K1 = spherogram.Link('10_74')
    
    # Compute the HOMFLY polynomial. In spherogram, the variables are l and m,
    # corresponding to the standard v and z respectively.
    P1_poly = K1.homfly_polynomial()
    
    # We need the v-span (l-span) of the polynomial.
    # We treat the polynomial as a polynomial in the variable 'l'.
    l = P1_poly.variables[0]
    
    # To find the lowest degree, we can substitute l with 1/l and find the new highest degree.
    P1_expanded = sympy.expand(P1_poly)
    P1_inv_l = P1_expanded.subs(l, 1/l)
    
    max_l_power = sympy.degree(P1_expanded, gen=l)
    min_l_power = -sympy.degree(sympy.expand(P1_inv_l), gen=l)
    
    v_span_K1 = max_l_power - min_l_power
    
    # Calculate the lower bound for the braid index (and min Seifert circles)
    # using the Morton-Franks-Williams inequality.
    lower_bound_K1 = 0.5 * (v_span_K1 + 1)

    # --- Part 2: Analyze K2 = closure of (sigma_1^-1)^3 * sigma_2^-1 ---
    
    # Define the braid word for K2. In spherogram, sigma_i is represented by i,
    # and its inverse by -i. The generators are 1-indexed.
    braid_word_K2 = [-1, -1, -1, -2]
    K2_link = spherogram.Link(braid=braid_word_K2)
    
    # Identify the knot. The identify() method returns a list of matching knots.
    # The result L4a1(0,0) is another name for the figure-eight knot, 4_1.
    identified_knot_name = K2_link.identify()[0].name()

    # Get the braid index for the identified knot.
    # The braid index of the figure-eight knot (4_1) is 3.
    braid_index_K2 = spherogram.Link(identified_knot_name).braid_index
    
    # --- Part 3: Calculate the difference ---
    
    difference = braid_index_K2 - lower_bound_K1
    
    # Print the final result in the required format.
    print(f"The braid index of K2 is: {braid_index_K2}")
    print(f"The lower bound for the minimum number of Seifert circles of K1 is: {lower_bound_K1}")
    print("The difference is calculated as (braid index of K2) - (lower bound for K1):")
    print(f"{braid_index_K2} - {lower_bound_K1} = {difference}")

solve_knot_problem()