import math

def solve_circulons():
    """
    Calculates the number of circulons for G=SO(3) in d+1 dimensions, for d=1..6.

    The number of circulons (circular defects in R^d) is given by the number of
    homotopy classes of maps from the complement of the circle, R^d \ S^1, to the
    gauge group G=SO(3).
    """

    # We need the orders of the homotopy groups of SO(3).
    # pi_k(SO(3)) = pi_k(S^3) for k>=2. pi_1(SO(3)) = Z_2.
    # We represent infinite order with math.inf. The trivial group has order 1.
    pi_so3_orders = {
        1: 2,          # pi_1(SO(3)) = Z_2
        2: 1,          # pi_2(SO(3)) = pi_2(S^3) = 0
        3: math.inf,   # pi_3(SO(3)) = pi_3(S^3) = Z
        4: 2,          # pi_4(SO(3)) = pi_4(S^3) = Z_2
        5: 2,          # pi_5(SO(3)) = pi_5(S^3) = Z_2
    }

    print("Calculating the number of circulons for G=SO(3) in d+1 dimensions.\n")

    # Case d=1
    d = 1
    # A circle (S^1) cannot be embedded in a line (R^1).
    # Thus, no such defects can exist.
    num_circulons_d1 = 0
    print(f"For d={d}: A circle cannot be embedded in R^1. The number of circulons is {num_circulons_d1}.")
    
    # Case d=2
    d = 2
    # The complement space R^2 \ S^1 has the homotopy type of S^1 v S^1.
    # The number of maps is |pi_1(SO(3))| * |pi_1(SO(3))|.
    order_pi1 = pi_so3_orders[1]
    num_circulons_d2 = order_pi1 * order_pi1
    print(f"For d={d}: The number is |pi_1(SO(3))| * |pi_1(SO(3))| = {order_pi1} * {order_pi1} = {num_circulons_d2}.")

    # Case d=3
    d = 3
    # The complement of a knot R^3 \ K is a K(G_K, 1) space.
    # Maps are classified by Hom(pi_1(R^3 \ K), pi_1(SO(3))) = Hom(G_K, Z_2).
    # The abelianization of any knot group G_K is Z, so |Hom(Z, Z_2)| = 2.
    num_circulons_d3 = 2
    print(f"For d={d}: The number is |Hom(pi_1(R^3 \\ K), pi_1(SO(3)))| = |Hom(G_K, Z_2)| = {num_circulons_d3}.")

    # Cases d > 3
    for d in range(4, 7):
        # The complement space R^d \ S^1 has the homotopy type of S^{d-1} v S^1.
        # The number of maps is |pi_{d-1}(SO(3))| * |pi_1(SO(3))|.
        k = d - 1
        order_pi_k = pi_so3_orders[k]
        order_pi1 = pi_so3_orders[1]
        
        if order_pi_k is math.inf:
            num_circulons = "infinite"
            print(f"For d={d}: The number is |pi_{k}(SO(3))| * |pi_1(SO(3))| = |Z| * |Z_2| = {num_circulons}.")
        else:
            num_circulons = order_pi_k * order_pi1
            print(f"For d={d}: The number is |pi_{k}(SO(3))| * |pi_1(SO(3))| = {order_pi_k} * {order_pi1} = {num_circulons}.")

solve_circulons()
<<<[0, 4, 2, "infinite", 4, 4]>>>