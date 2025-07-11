import sympy

def find_possible_degrees():
    """
    This function calculates the possible degrees of normal field extensions L/Q
    for the splitting field K of the given polynomial, where Q < L < K.
    """
    x = sympy.Symbol('x')
    P = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    # Step 1: Factor the polynomial over the rational numbers.
    # The result shows P(x) = (x**2 - 2) * (x**5 - 9*x + 3).
    quad_factor = x**2 - 2
    quintic_factor = x**5 - 9*x + 3
    
    # Step 2 & 3: Determine the Galois group G = Gal(K/Q).
    # The Galois group of the quadratic factor x**2 - 2 is the cyclic group C_2.
    # The quintic factor x**5 - 9*x + 3 is irreducible over Q (by Eisenstein's
    # criterion with p=3). It has exactly 3 real roots, so its Galois group is S_5.
    # To confirm the overall Galois group is the direct product C_2 x S_5, we must
    # ensure the splitting fields K1=Q(sqrt(2)) and K2 (of the quintic) only
    # intersect in Q. This is true if K2 does not contain sqrt(2). The only quadratic
    # subfield of K2 is Q(sqrt(D)) where D is the discriminant of the quintic.
    # D = 14863419. Since D is not of the form 2*k^2, Q(sqrt(D)) != Q(sqrt(2)).
    # Thus, G = C_2 x S_5.
    
    # The order of G is |C_2| * |S_5|.
    order_C2 = 2
    order_S5 = sympy.factorial(5)
    order_G = order_C2 * order_S5 # 2 * 120 = 240

    # Step 4 & 5: Find proper non-trivial normal subgroups H of G and their indices.
    # The normal subgroups of S_5 are {e}, A_5, S_5.
    # The normal subgroups of C_2 are {e}, C_2.
    # This gives normal subgroups of G of the form H1 x H2.
    # We also have a "diagonal" normal subgroup because S_5/A_5 is isomorphic to C_2/{e}.
    
    order_A5 = order_S5 // 2 # 120 / 2 = 60

    normal_subgroup_orders = [
        order_C2 * 1,         # C_2 x {e}
        1 * order_A5,         # {e} x A_5
        1 * order_S5,         # {e} x S_5
        order_C2 * order_A5,  # C_2 x A_5
        order_S5              # The diagonal subgroup N = {(c,s) | c = sgn(s)} has order |S_5|=120
    ]
    
    # Calculate the corresponding degrees [L:Q] = |G|/|H|.
    possible_degrees = set()
    for h_order in normal_subgroup_orders:
        degree = order_G // h_order
        possible_degrees.add(degree)

    # Step 6: Print the sorted list of unique degrees.
    final_degrees = sorted(list(possible_degrees))
    print("The possible degrees of the normal field extensions are:")
    # We unpack the list to print the numbers.
    print(*final_degrees, sep=', ')

find_possible_degrees()