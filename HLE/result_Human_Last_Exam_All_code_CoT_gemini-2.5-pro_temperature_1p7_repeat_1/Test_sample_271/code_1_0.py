import math

def solve_degrees():
    """
    Calculates and lists the possible degrees of normal field extensions L/Q
    for the splitting field K of the polynomial f(x) over Q, where
    Q is a proper subfield of L and L is a proper subfield of K.
    """
    
    # Step 1: Factor the polynomial and identify the components.
    # The polynomial is f(x) = x^7 - 2x^5 - 9x^3 + 3x^2 + 18x - 6.
    # This polynomial factors over Q into:
    # f(x) = (x^2 - 2) * (x^5 - 9x + 3)
    # Let f1(x) = x^2 - 2 and f2(x) = x^5 - 9x + 3.

    # Step 2: Determine the Galois group for each factor.
    # The Galois group of f1(x) = x^2 - 2 over Q is the cyclic group C_2.
    g1_order = 2
    
    # The polynomial f2(x) = x^5 - 9x + 3 is irreducible over Q by Eisenstein's criterion for p=3.
    # It can be shown to have exactly 3 real roots and 2 complex conjugate roots.
    # A key theorem in Galois theory states that an irreducible polynomial of prime degree p with
    # exactly p-2 real roots has the symmetric group S_p as its Galois group. Here p=5.
    # So, the Galois group of f2(x) over Q is S_5.
    g2_order = math.factorial(5)  # |S_5| = 120

    # Step 3: Determine the full Galois group G.
    # The splitting field K is the compositum of the splitting fields of f1 and f2.
    # The splitting field of f1 is Q(sqrt(2)).
    # The splitting field of f2, K2, has a unique quadratic subfield Q(sqrt(D)), where D is the discriminant.
    # D = 256*(-9)**5 + 3125*(3)**4 = -14863419.
    # Since Q(sqrt(2)) is not isomorphic to Q(sqrt(-14863419)), the intersection of the splitting fields is Q.
    # Therefore, the full Galois group G is the direct product C_2 x S_5.
    G_order = g1_order * g2_order
    print(f"The Galois group is G = C2 x S5, and its order is {g1_order} * {g2_order} = {G_order}.")

    # Step 4: Identify the proper non-trivial normal subgroups H of G.
    # The proper non-trivial normal subgroups of G = C_2 x S_5 are:
    # 1. {e} x A5 (A5 is the alternating group, a normal subgroup of S5)
    # 2. {e} x S5
    # 3. C2 x {e}
    # 4. C2 x A5
    # 5. A "diagonal" subgroup, call it N_d, which exists because C2 is isomorphic to the quotient group S5/A5.
    #
    # The orders of these subgroups are:
    A5_order = g2_order // 2
    H_orders = {
        "{e} x A5": A5_order,
        "{e} x S5": g2_order,
        "C2 x {e}": g1_order,
        "C2 x A5": g1_order * A5_order,
        "N_d": g1_order * A5_order,
    }

    # Step 5: Calculate the index |G:H| for each subgroup. This gives the degree of the extension.
    print("\nThe degrees of the normal intermediate fields are the indices of the proper normal subgroups:")
    possible_degrees = set()
    for name, order in H_orders.items():
        index = G_order // order
        print(f"For subgroup H = {name} with order {order}, the degree is |G|/|H| = {G_order} / {order} = {index}")
        possible_degrees.add(index)
        
    # Step 6: List the unique degrees found.
    print("\nThe set of all possible degrees for these normal field extensions is:")
    sorted_degrees = sorted(list(possible_degrees))
    print(', '.join(map(str, sorted_degrees)))

solve_degrees()