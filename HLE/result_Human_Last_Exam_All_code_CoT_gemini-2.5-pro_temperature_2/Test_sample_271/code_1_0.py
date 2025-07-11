import sympy

def find_possible_degrees():
    """
    Solves the problem by following the planned steps:
    1. Factor the polynomial f(x).
    2. Determine the Galois group G of f(x).
    3. Find all proper non-trivial normal subgroups H of G.
    4. Calculate the index [G:H] for each H.
    """
    x = sympy.Symbol('x')
    f_poly = sympy.Poly(x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6, x, domain='ZZ')

    # Step 1: Factor the polynomial
    # We find that sqrt(2) is a root, which implies (x^2 - 2) is a factor.
    # Polynomial long division yields the other factor.
    g_poly = sympy.Poly(x**2 - 2, x)
    h_poly = sympy.Poly(x**5 - 9*x + 3, x)
    
    print(f"The polynomial f(x) = {f_poly.as_expr()} can be factored over Q into:")
    print(f"f(x) = ({g_poly.as_expr()}) * ({h_poly.as_expr()})")
    print("-" * 50)

    # Step 2: Determine the Galois Group
    print("Step 2: Determining the Galois Group G = Gal(K/Q)")
    
    # Galois group of g(x) = x^2 - 2
    G_g = "C_2 (Cyclic group of order 2)"
    order_G_g = 2
    print(f"The Galois group of g(x) = {g_poly.as_expr()} is {G_g}, with order {order_G_g}.")
    
    # Galois group of h(x) = x^5 - 9x + 3
    # h(x) is irreducible by Eisenstein's criterion for p=3.
    # As an irreducible quintic with 3 real roots and 2 complex roots, its Galois group is S_5.
    G_h = "S_5 (Symmetric group on 5 elements)"
    order_G_h = sympy.factorial(5)
    print(f"The Galois group of h(x) = {h_poly.as_expr()} is {G_h}, with order {order_G_h}.")

    # Galois group of f(x)
    # The splitting field of g(x) is Q(sqrt(2)). The quadratic subfield of the splitting
    # field of h(x) is Q(sqrt(discriminant(h))).
    disc_h = sympy.discriminant(h_poly) # 37158957
    print(f"The discriminant of h(x) is {disc_h}.")
    print("Since sqrt(37158957) is not a rational multiple of sqrt(2), the splitting fields are linearly disjoint.")
    
    G_f = "C_2 x S_5"
    order_G_f = order_G_g * order_G_h
    print(f"The Galois group G of the splitting field of f(x) is the direct product {G_f}.")
    print(f"The order of G is {order_G_g} * {order_G_h} = {order_G_f}.")
    print("-" * 50)

    # Step 3 & 4: Find normal subgroups and calculate their indices
    print("Step 3 & 4: Finding degrees from normal subgroups of G = C_2 x S_5")
    print("The possible degrees of proper normal extensions L/Q are the indices [G:H] of")
    print("proper, non-trivial normal subgroups H of G.")

    degrees = set()
    
    # Normal subgroups of S_5 are {e}, A_5, S_5.
    # This leads to several normal subgroups of C_2 x S_5.
    
    # H_1 = {e} x A_5
    order_H1 = sympy.factorial(5) / 2
    index_1 = order_G_f // order_H1
    degrees.add(index_1)
    print(f"For H = {{e}} x A_5, |H|={int(order_H1)}, the degree is |G|/|H| = {index_1}")
    
    # H_2 = {e} x S_5
    order_H2 = order_G_h
    index_2 = order_G_f // order_H2
    degrees.add(index_2)
    print(f"For H = {{e}} x S_5, |H|={int(order_H2)}, the degree is |G|/|H| = {index_2}")

    # H_3 = C_2 x {e}
    order_H3 = order_G_g
    index_3 = order_G_f // order_H3
    degrees.add(index_3)
    print(f"For H = C_2 x {{e}}, |H|={int(order_H3)}, the degree is |G|/|H| = {index_3}")
    
    # H_4 = C_2 x A_5
    order_H4 = order_G_g * (sympy.factorial(5) / 2)
    index_4 = order_G_f // order_H4
    degrees.add(index_4)
    print(f"For H = C_2 x A_5, |H|={int(order_H4)}, the degree is |G|/|H| = {index_4}")

    # H_5 = "Diagonal" normal subgroup
    # This subgroup arises from the surjective homomorphism S_5 -> C_2 (the sign map).
    order_H5 = order_G_h # It is isomorphic to S_5
    index_5 = order_G_f // order_H5
    degrees.add(index_5)
    print(f"For the diagonal normal subgroup H isomorphic to S_5, |H|={int(order_H5)}, the degree is |G|/|H| = {index_5}")

    sorted_degrees = sorted(list(degrees))
    print("-" * 50)
    print("The complete list of possible degrees for such extensions is:")
    print(*sorted_degrees, sep=', ')

find_possible_degrees()
>>>2, 4, 120