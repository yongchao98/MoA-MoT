import math

def solve():
    """
    Calculates the maximal chromatic number of a graph G that is the sum of three cycles of length n.
    """
    # The chromatic number of a graph G = G1 + G2 + ... + Gk is chi(G) = chi(G1) * chi(G2) * ... * chi(Gk).
    # In our case, G = C_n + C_n + C_n, so chi(G) = chi(C_n)^3.

    # The chromatic number of a cycle C_n depends on n:
    # chi(C_n) = 2 if n is even.
    # chi(C_n) = 3 if n is odd.

    # To find the maximal chromatic number of G, we must use the maximal
    # possible value for chi(C_n).
    chi_cn_even = 2
    chi_cn_odd = 3
    max_chi_cn = max(chi_cn_even, chi_cn_odd)

    # The number of cycles being summed.
    num_cycles = 3

    # The maximal chromatic number of G is (max_chi_cn) ^ num_cycles.
    maximal_chi_g = math.pow(max_chi_cn, num_cycles)

    print("The maximal chromatic number of G is found by maximizing chi(C_n).")
    print(f"The maximum chromatic number for a single cycle C_n is {max_chi_cn} (when n is odd).")
    print(f"The graph G is the sum of {num_cycles} such cycles.")
    print("\nThe final equation for the maximal chromatic number of G is:")
    # Print each number in the final equation as requested.
    print(f"{max_chi_cn} ^ {num_cycles} = {int(maximal_chi_g)}")

solve()