import math

def euler_phi(n):
    """Computes Euler's totient function phi(n)."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def get_divisors(n):
    """Computes all divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve():
    """
    Calculates the number of non-isomorphic, connected, 3-regular adjustable graphs
    with 2000 vertices and a perfect matching.
    """
    n = 1000

    # The number of such graphs is equivalent to the number of binary bracelets of length n.
    # The formula is given by Burnside's Lemma applied to the dihedral group D_n acting on
    # the vertices of a regular n-gon, with 2 colors.
    # N = (1 / 2n) * (sum_{d|n} phi(d)*2^(n/d) + reflection_term)
    
    # Term from rotations in the dihedral group
    divs = get_divisors(n)
    rotation_term = 0
    for d in divs:
        term = euler_phi(d) * (1 << (n // d))
        rotation_term += term
    
    # Term from reflections in the dihedral group (for even n)
    # The formula is (n/2)*2^((n+2)/2) + (n/2)*2^(n/2)
    term1_ref = (n // 2) * (1 << ((n + 2) // 2))
    term2_ref = (n // 2) * (1 << (n // 2))
    reflection_term = term1_ref + term2_ref

    # Total sum in the numerator of Burnside's Lemma
    total_sum = rotation_term + reflection_term
    
    # The order of the dihedral group D_n is 2n
    group_order = 2 * n
    
    # The result must be an integer
    num_graphs = total_sum // group_order
    
    print(f"The number of graphs is calculated using the formula for counting binary bracelets of length {n}.")
    print("Formula: N = (1 / (2*n)) * [ (sum over d|n of phi(d)*2^(n/d)) + reflection_term ]")
    print(f"\nFor n = {n}:")
    print(f"The sum from rotations is: {rotation_term}")
    print(f"The sum from reflections is: {reflection_term}")
    print(f"The total numerator is: {total_sum}")
    print(f"The group order (denominator) is: {group_order}")
    print(f"\nThe total number of non-isomorphic graphs is: {num_graphs}")

solve()