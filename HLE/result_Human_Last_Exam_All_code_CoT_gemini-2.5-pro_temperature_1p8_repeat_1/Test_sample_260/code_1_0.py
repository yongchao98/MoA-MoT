import math

def solve():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Step 1: Define the properties of the base spaces X_i.
    
    # X1 is from a pentagon (5 sides). pi_1(X1) is Z_5.
    # Its CW-structure is one 0-cell, one 1-cell, one 2-cell.
    order_pi1_X1 = 5
    chi_X1 = 1 - 1 + 1

    # X2 is from an octagon (8 sides). pi_1(X2) is Z_8.
    # Its CW-structure is one 0-cell, one 1-cell, one 2-cell.
    order_pi1_X2 = 8
    chi_X2 = 1 - 1 + 1
    
    # X3 is the real projective plane. pi_1(X3) is Z_2.
    # Its standard CW-structure is one 0-cell, one 1-cell, one 2-cell.
    order_pi1_X3 = 2
    chi_X3 = 1 - 1 + 1
    
    orders = [order_pi1_X1, order_pi1_X2, order_pi1_X3]
    chis = [chi_X1, chi_X2, chi_X3]
    n = len(orders)
    
    # Step 2: Calculate properties of the connected sum Y.
    
    # The first homology group H_1(Y) is the direct sum H_1(X1) + H_1(X2) + H_1(X3).
    # H_1(X_i) is the abelianization of pi_1(X_i), which is Z_n for pi_1(X_i) = Z_n.
    # So H_1(Y) is isomorphic to Z_5 + Z_8 + Z_2.
    # The order of H_1(Y) is the product of the orders of the component groups.
    order_H1_Y = math.prod(orders)
    
    # The Euler characteristic of the connected sum Y = X1 # X2 # X3 is
    # chi(Y) = chi(X1) + chi(X2) + chi(X3) - 2*(n-1).
    chi_Y = sum(chis) - 2 * (n - 1)
    
    # Step 3: Use covering space theory to find the rank.
    
    # The kernel K = Ker(h_*) is the commutator subgroup of pi_1(Y).
    # The degree of the covering space Y_K corresponding to K is |H_1(Y)|.
    # The Euler characteristic of the covering space is chi(Y_K) = |H_1(Y)| * chi(Y).
    chi_YK = order_H1_Y * chi_Y
    
    # K = pi_1(Y_K) is a free group of rank r. Y_K is homotopy equivalent
    # to a wedge of r circles, whose Euler characteristic is 1-r.
    # So, chi(Y_K) = 1 - r, which means r = 1 - chi(Y_K).
    rank_K = 1 - chi_YK
    
    # Step 4: Print the results.
    print("This script calculates the rank of the kernel of the Hurewicz homomorphism h*: pi_1(Y) -> H_1(Y).")
    print(f"The space Y is the connected sum of X1, X2, and X3.")
    print("-" * 30)
    
    print(f"1. Order of the first homology group H_1(Y):")
    print(f"|H_1(Y)| = |H_1(X1)| * |H_1(X2)| * |H_1(X3)| = {order_pi1_X1} * {order_pi1_X2} * {order_pi1_X3} = {order_H1_Y}")
    
    print("\n2. Euler characteristic of Y:")
    print(f"chi(Y) = chi(X1) + chi(X2) + chi(X3) - 2*(n-1) = {chi_X1} + {chi_X2} + {chi_X3} - 2*({n}-1) = {chi_Y}")

    print("\n3. Rank of the kernel K:")
    print("The rank 'r' is calculated using the formula r = 1 - |H_1(Y)| * chi(Y).")
    print("\nThe final equation is:")
    print(f"rank = 1 - (|H_1(Y)| * chi(Y))")
    print(f"rank = 1 - ({order_H1_Y} * {chi_Y})")
    print(f"rank = 1 - ({chi_YK})")
    print(f"rank = {rank_K}")

solve()