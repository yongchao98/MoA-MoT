import math

def count_power_subgroups():
    """
    Calculates the number of power subgroups in the semidihedral group of size 512.

    A power subgroup of a group G is a subgroup H that can be expressed as
    H = G^m = {g^m | g in G} for some integer m.

    For the semidihedral group G = SD_{2^n}, G^m is a subgroup if and only if:
    1. m is even.
    2. m is congruent to 1 modulo 2^(n-2).

    We analyze these two cases for G = SD_512, where n=9.
    """
    n = 9
    order_r = 2**(n - 1)  # This is 256, the order of the maximal cyclic subgroup <r>

    # Case 1: m is an even integer.
    # The resulting power subgroups are the cyclic subgroups <r^d> where d is an even divisor of 256.
    # The divisors of 256 (which is 2^8) are 1, 2, 4, 8, 16, 32, 64, 128, 256.
    # The even divisors are 2, 4, 8, 16, 32, 64, 128, 256. There are 8 of them.
    # This number can be calculated as log2 of the order of <r>.
    num_even_case_subgroups = int(math.log2(order_r))

    # Case 2: m is congruent to 1 modulo 2^(n-2).
    # For n=9, this condition is m = 1 (mod 128).
    # This always results in the entire group G = SD_512 itself.
    # This gives exactly one subgroup.
    num_odd_case_subgroups = 1

    # The subgroups from Case 1 are all cyclic and proper subgroups of G.
    # The subgroup from Case 2 is G itself, which is non-cyclic.
    # Therefore, the two sets of subgroups are disjoint.
    total_subgroups = num_even_case_subgroups + num_odd_case_subgroups

    print(f"The number of power subgroups from Case 1 (m is even) is {num_even_case_subgroups}.")
    print(f"The number of power subgroups from Case 2 (m = 1 mod 128) is {num_odd_case_subgroups}.")
    print("The total number of power subgroups is the sum of these two disjoint sets.")
    print("Final Calculation:")
    print(f"{num_even_case_subgroups} + {num_odd_case_subgroups} = {total_subgroups}")


count_power_subgroups()