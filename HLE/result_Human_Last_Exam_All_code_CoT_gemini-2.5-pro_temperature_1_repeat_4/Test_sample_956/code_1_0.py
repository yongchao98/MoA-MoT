import math
from collections import Counter
from functools import reduce
import operator
from sympy.combinatorics import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup

def count_iso_p_subgroups(p, lam):
    """
    Counts the number of non-isomorphic subgroups of an abelian p-group
    of type lambda. Lambda is a partition (a list of integer exponents).
    This function implements known values for small, common partitions.
    """
    lam_tuple = tuple(sorted(lam, reverse=True))
    if not lam_tuple: # Trivial group component
        return 1
    
    type_str = " x ".join([f"C_{p**k}" for k in lam_tuple])
    print(f"The Sylow {p}-subgroup is isomorphic to {type_str}, with partition type {lam_tuple}.")
    
    if len(lam_tuple) == 1: # Cyclic group C_{p^k}
        k = lam_tuple[0]
        count = k + 1
        print(f"For a cyclic p-group, the number of non-isomorphic subgroups is k+1 = {k} + 1 = {count}")
        return count
    if lam_tuple == (1, 1): # C_p x C_p
        count = 3
        print(f"For a group of type (1, 1), the number of non-isomorphic subgroups is {count}.")
        return count
    if lam_tuple == (2, 1): # C_{p^2} x C_p
        count = 5
        print(f"For a group of type (2, 1), the number of non-isomorphic subgroups is {count}.")
        return count
    if lam_tuple == (1, 1, 1): # C_p x C_p x C_p
        count = 5
        print(f"For a group of type (1, 1, 1), the number of non-isomorphic subgroups is {count}.")
        return count

    # Fallback for unimplemented partitions
    print(f"Warning: Calculation for partition type {lam_tuple} is not implemented.")
    return -1 

def get_prime_factorization(n):
    """Computes the prime factorization of n as a Counter."""
    factors = Counter()
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def solve():
    """
    Main function to solve the problem.
    """
    # The set is {1, 2, 3, 4, 5, 6, 7, 8, 9, x, y, z}.
    # We map it to {0, 1, ..., 11} for sympy.
    # 1->0, 2->1, 3->2, 4->3, 5->4, 6->5, 7->6, 8->7, 9->8, x->9, y->10, z->11
    
    # a = (1, 3, 2, 5, 4, 7, 6, 9, 8, y, x) -> (0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9)
    # b = (1, 8, 5, 9)(4, x, 7, 6) -> (0, 7, 4, 8)(3, 9, 6, 5)
    # c = (1, 2)(3, z)(4, 8)(5, 6)(7, y)(9, x) -> (0, 1)(2, 11)(3, 7)(4, 5)(6, 10)(8, 9)

    a_perm = Permutation([0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9])
    b_perm = Permutation([[0, 7, 4, 8], [3, 9, 6, 5]])
    c_perm = Permutation([[0, 1], [2, 11], [3, 7], [4, 5], [6, 10], [8, 9]])

    G = PermutationGroup([a_perm, b_perm, c_perm])

    A_factors = G.schur_multiplier()

    if not A_factors or (len(A_factors) == 1 and A_factors[0] == 1):
        print("The Schur multiplier A is the trivial group, A = C_1.")
        print("The only subgroup is C_1 itself, which is not proper.")
        print("The number of proper subgroups of A up to isomorphism is 0.")
        return

    A_structure_str = " x ".join([f"C_{n}" for n in A_factors])
    print(f"The Schur multiplier A is isomorphic to {A_structure_str}.")

    all_prime_factors = set()
    for n in A_factors:
        for p in get_prime_factorization(n):
            all_prime_factors.add(p)

    num_iso_subgroups_per_p = []
    num_iso_subgroups_per_p_str = []

    for p in sorted(list(all_prime_factors)):
        lambda_p = []
        for n in A_factors:
            exponent = get_prime_factorization(n)[p]
            if exponent > 0:
                lambda_p.append(exponent)
        
        num_p_subgroups = count_iso_p_subgroups(p, lambda_p)
        if num_p_subgroups == -1: return

        num_iso_subgroups_per_p.append(num_p_subgroups)
        num_iso_subgroups_per_p_str.append(str(num_p_subgroups))

    if len(num_iso_subgroups_per_p) > 1:
        total_iso_subgroups = reduce(operator.mul, num_iso_subgroups_per_p)
        print(f"\nTotal number of non-isomorphic subgroups of A is the product of these counts: "
              f"{' * '.join(num_iso_subgroups_per_p_str)} = {total_iso_subgroups}")
    else:
        total_iso_subgroups = num_iso_subgroups_per_p[0]
        print(f"\nTotal number of non-isomorphic subgroups of A is {total_iso_subgroups}.")

    num_proper_iso_subgroups = total_iso_subgroups - 1
    print(f"The number of proper subgroups of A up to isomorphism is the total minus 1 (for A itself): "
          f"{total_iso_subgroups} - 1 = {num_proper_iso_subgroups}")
    print(f"\n<<< {num_proper_iso_subgroups} >>>")

if __name__ == '__main__':
    solve()