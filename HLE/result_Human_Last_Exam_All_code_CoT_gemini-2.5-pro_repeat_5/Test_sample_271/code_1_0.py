from sympy.combinatorics.perm_groups import SymmetricGroup, CyclicGroup
from sympy.combinatorics.group_constructs import DirectProduct

def solve_degrees():
    """
    This function calculates the possible degrees of intermediate normal field extensions.
    """
    # 1. Define the Galois group G = S5 x C2.
    S5 = SymmetricGroup(5)
    C2 = CyclicGroup(2)
    G = DirectProduct(S5, C2)
    order_G = G.order()

    print(f"The Galois group G of the polynomial is the direct product S5 x C2.")
    print(f"The order of G is {order_G}.")
    print("\nThe normal extensions L/Q with Q < L < K correspond to proper non-trivial normal subgroups H of G.")
    print("The degree of such an extension is [L:Q] = |G|/|H|.\n")
    
    # 2. Find all normal subgroups of G.
    normal_subgroups = G.normal_subgroups()
    
    # 3. Calculate the indices for proper, non-trivial normal subgroups.
    possible_degrees = set()
    
    print("Calculating the degrees based on the orders of normal subgroups:")
    for H in normal_subgroups:
        order_H = H.order()
        # We are interested in proper, non-trivial normal subgroups.
        if 1 < order_H < order_G:
            degree = order_G // order_H
            possible_degrees.add(degree)
            print(f"Found a proper normal subgroup H with order {order_H}. The degree of the extension is {order_G} / {order_H} = {degree}.")

    # 4. List all the unique possible degrees.
    sorted_degrees = sorted(list(possible_degrees))
    print(f"\nThe list of all possible degrees for the normal field extensions L/Q is: {sorted_degrees}")

solve_degrees()