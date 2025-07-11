import math

def solve_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for Y = X1 # X2 # X3.
    """
    # Step 1: Define the orders of the fundamental groups of X1, X2, and X3.
    # pi_1(X1) = Z_5, pi_1(X2) = Z_8, pi_1(X3) = Z_2
    n1 = 5
    n2 = 8
    n3 = 2
    orders = [n1, n2, n3]
    num_groups = len(orders)

    print("Step 1: The fundamental groups are G1 = Z_5, G2 = Z_8, and G3 = Z_2.")
    print(f"Their orders are |G1| = {n1}, |G2| = {n2}, |G3| = {n3}.\n")

    # Step 2: Calculate the Euler characteristic of the fundamental group G = G1 * G2 * G3.
    # chi(G) = sum(chi(Gi)) - (num_groups - 1)
    # chi(Gi) = 1/|Gi|
    chi_G = sum(1/n for n in orders) - (num_groups - 1)
    
    print("Step 2: The fundamental group of the connected sum is G = G1 * G2 * G3.")
    print("The Euler characteristic of G is calculated as:")
    print(f"chi(G) = (1/{n1} + 1/{n2} + 1/{n3}) - ({num_groups} - 1)")
    print(f"chi(G) = {sum(1/n for n in orders)} - {num_groups - 1} = {chi_G}\n")

    # Step 3: Calculate the Euler characteristic of the first homology group H1(Y).
    # H1(Y) = G_ab = G1 + G2 + G3 (direct sum)
    # |H1(Y)| = |G1| * |G2| * |G3|
    order_H1 = math.prod(orders)
    chi_H1 = 1 / order_H1

    print("Step 3: The first homology group H1(Y) is the abelianization of G, which is H1 = G1 + G2 + G3.")
    print(f"The order of H1 is |H1| = {n1} * {n2} * {n3} = {order_H1}.")
    print(f"The Euler characteristic of H1 is chi(H1) = 1 / |H1| = 1/{order_H1} = {chi_H1}\n")

    # Step 4: Calculate the rank of the kernel K.
    # The kernel K is a free group of rank r, so chi(K) = 1 - r.
    # From the short exact sequence 1 -> K -> G -> H1 -> 1, we have chi(G) = chi(K) * chi(H1).
    # So, chi(K) = chi(G) / chi(H1).
    # 1 - r = chi(G) / chi(H1)
    # r = 1 - (chi(G) / chi(H1))
    one_minus_r = chi_G / chi_H1
    rank_r = 1 - one_minus_r

    print("Step 4: The kernel K is a free group of rank r. Its Euler characteristic is chi(K) = 1 - r.")
    print("Using the formula chi(G) = chi(K) * chi(H1), we solve for r.")
    print(f"1 - r = chi(G) / chi(H1)")
    print(f"1 - r = {chi_G} / {chi_H1} = {one_minus_r}")
    print(f"r = 1 - ({one_minus_r})")
    print(f"The final rank is r = {int(rank_r)}.")

solve_rank()