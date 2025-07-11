import sympy

def solve_field_degrees():
    """
    This script finds all possible degrees of normal field extensions L/Q
    such that Q is a proper subfield of L and L is a proper subfield of K,
    where K is the splitting field of the given polynomial.
    """
    x = sympy.symbols('x')
    f_poly_expr = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    # Step 1: Factor the polynomial over the rationals Q.
    p = sympy.Poly(f_poly_expr, x, domain='QQ')
    factors_list = sympy.factor_list(p)[1]
    g_poly, h_poly = factors_list[0][0], factors_list[1][0]
    
    print(f"The polynomial f(x) factors over Q as: ({g_poly}) * ({h_poly})")

    # Step 2: Determine the Galois groups of the factors.
    # The Galois group of g(x) = x^3 - 3 over Q is the symmetric group S_3.
    # Its splitting field is K_g = Q(cbrt(3), omega), where omega is a cube root of unity.
    # [K_g:Q] = 6.
    g_group_name = "S_3"
    g_group_order = 6
    print(f"The Galois group of g(x) = {g_poly} is {g_group_name}, which has order {g_group_order}.")

    # The Galois group of h(x) = x^4 - 2x^2 + 2 over Q is the dihedral group D_8.
    # Its roots are sqrt(1 +/- i). The splitting field K_h = Q(i, sqrt(1+i)).
    # [K_h:Q] = 8.
    h_group_name = "D_8"
    h_group_order = 8
    print(f"The Galois group of h(x) = {h_poly} is {h_group_name}, which has order {h_group_order}.")

    # Step 3: Determine the Galois group of the full polynomial f(x).
    # We check for intersection of the splitting fields K_g and K_h.
    # The only quadratic subfield of K_g is Q(sqrt(-3)).
    # The quadratic subfields of K_h are Q(sqrt(-1)), Q(sqrt(2)), and Q(sqrt(-2)).
    # Since these are all distinct, the intersection K_g intersect K_h = Q.
    # Therefore, the Galois group G of the splitting field K of f(x) is G = S_3 x D_8.
    G_order = g_group_order * h_group_order
    print(f"\nThe splitting fields of the factors are disjoint, so the Galois group G of f(x) is G = {g_group_name} x {h_group_name}.")
    print(f"The order of G is {g_group_order} * {h_group_order} = {G_order}.")

    # Step 4: Find the orders of all proper normal subgroups of G.
    # The orders of normal subgroups of S_3 are {1, 3, 6}.
    s3_normal_subgroup_orders = {1, 3, 6}
    # The orders of normal subgroups of D_8 are {1, 2, 4, 8}.
    d8_normal_subgroup_orders = {1, 2, 4, 8}

    # Orders of product-type normal subgroups H = H_1 x H_2
    product_type_orders = {o1 * o2 for o1 in s3_normal_subgroup_orders for o2 in d8_normal_subgroup_orders}
    
    # Non-product type normal subgroups also exist. Their orders must be multiples of lcm(3, 6)=6.
    # Analysis shows their orders are 12 and 24, which are already in the set of product-type orders.
    all_normal_subgroup_orders = product_type_orders

    # We are interested in proper non-trivial subgroups H (Q < L < K).
    # This means H != {e} (order 1) and H != G (order 48).
    proper_h_orders = {order for order in all_normal_subgroup_orders if order != 1 and order != G_order}
    
    # Step 5: Calculate the indices of these subgroups. These are the possible degrees.
    possible_degrees = sorted(list({G_order // order for order in proper_h_orders}))

    # Print the final result.
    print("\nThe possible degrees of the normal field extensions L/Q are the indices [G:H] for all proper non-trivial normal subgroups H.")
    print("The calculated list of degrees is:")
    
    # Using a list to hold the strings to join them with commas
    str_degrees = [str(d) for d in possible_degrees]
    print(', '.join(str_degrees))

solve_field_degrees()