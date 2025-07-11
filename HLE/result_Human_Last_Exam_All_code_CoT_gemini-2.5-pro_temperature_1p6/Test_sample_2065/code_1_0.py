import math

def solve_circulons():
    """
    Calculates the number of circulon defects for G=SO(3) in d=1 to 6 spatial dimensions.
    """

    # We use pi_k(SO(3)). For k>=2, pi_k(SO(3)) = pi_k(S^3).
    # Using math.inf for Z, the group of integers.
    pi = {
        0: 1,           # pi_0(SO(3)) = 1 (path-connected)
        1: 2,           # pi_1(SO(3)) = Z_2
        2: 1,           # pi_2(SO(3)) = pi_2(S^3) = 0
        3: math.inf,    # pi_3(SO(3)) = pi_3(S^3) = Z
        4: 2,           # pi_4(SO(3)) = pi_4(S^3) = Z_2
        5: 2,           # pi_5(SO(3)) = pi_5(S^3) = Z_2
        6: 12           # pi_6(SO(3)) = pi_6(S^3) = Z_12
    }

    print("The number of circulons for gauge group G=SO(3) in d spatial dimensions:")

    for d in range(1, 7):
        if d == 1:
            # A circle cannot be embedded in R^1. The classification is trivial.
            result = 1
            print(f"For d=1: The defect cannot be embedded. The number of configurations is {result}.")
        
        elif d == 2:
            # The classifying space is S^1 u S^1.
            # N = |[S^1, SO(3)] x [S^1, SO(3)]| = |pi_1(SO(3))|^2
            pi_1 = pi[1]
            result = pi_1 ** 2
            print(f"For d=2: The number is |pi_1(SO(3))|^2 = {pi_1}^2 = {result}.")

        elif d == 3:
            # The classifying space is the torus T^2 = S^1 x S^1.
            # N = |[T^2, SO(3)]| = |H^1(T^2, Z_2)| = |Z_2 x Z_2|
            result = 4
            print(f"For d=3: The number is |H^1(T^2, Z_2)| = |Z_2 x Z_2| = {result}.")

        elif d > 3:
            # The classifying space is S^1 x S^{d-2}.
            # The number of circulons N(d) = |H^1(S^1 x S^{d-2}; Z_2)| * |pi_{d-2}(SO(3))| * |pi_{d-1}(SO(3))|
            # Note that for d > 3, d-2 >= 2, so pi_{d-2} and pi_{d-1} of SO(3) and S^3 are the same.
            # The H^1 term gives a factor of 2.
            factor_H1 = 2
            pi_d_minus_2 = pi.get(d - 2, 0)
            pi_d_minus_1 = pi.get(d - 1, 0)
            result = factor_H1 * pi_d_minus_2 * pi_d_minus_1
            
            # Format infinity for printing
            pi_d2_str = 'infinity' if pi_d_minus_2 == math.inf else pi_d_minus_2
            pi_d1_str = 'infinity' if pi_d_minus_1 == math.inf else pi_d_minus_1
            res_str = 'infinity' if result == math.inf else int(result)

            print(f"For d={d}: The number is 2 * |pi_{d-2}(SO(3))| * |pi_{d-1}(SO(3))| = {factor_H1} * {pi_d2_str} * {pi_d1_str} = {res_str}.")

solve_circulons()