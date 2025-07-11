import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    """
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
    """
    Gets all divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve():
    """
    Calculates the number of non-isomorphic graphs based on the described model.
    """
    n = 1000
    num_colors = 2
    group_order = 2 * n

    # --- Sum over rotations using Burnside's Lemma ---
    # The number of colorings fixed by a rotation of k steps is num_colors^gcd(k, n).
    # The sum over all rotations is sum_{d|n} phi(d) * (num_colors)^(n/d).
    # This is equivalent to sum_{d|n} phi(n/d) * (num_colors)^d.
    divisors = get_divisors(n)
    rotation_sum = 0
    
    print("The number of graphs is found by counting 2-colorings of a cycle C_1000 under the dihedral group D_1000.")
    print("Using Burnside's Lemma, N = (1/|G|) * sum_{g in G} |X^g|.")
    print(f"Here, |G| = 2*n = {group_order}, and the set X is the {num_colors**n} possible colorings.")
    
    print("\n1. Contribution from rotations:")
    print(f"Sum = sum over divisors d of n: phi(n/d) * {num_colors}^d")
    
    for d in divisors:
        term = phi(n // d) * (num_colors**d)
        rotation_sum += term
    
    print(f"Sum for all {n} rotations = {rotation_sum}")

    # --- Sum over reflections ---
    # For n even, there are n/2 reflections with num_colors^(n/2) fixed colorings,
    # and n/2 reflections with num_colors^(n/2 + 1) fixed colorings.
    num_refl_type1 = n // 2
    num_refl_type2 = n // 2
    
    fixed_type1 = num_colors**(n // 2)
    fixed_type2 = num_colors**(n // 2 + 1)
    
    reflection_sum = num_refl_type1 * fixed_type1 + num_refl_type2 * fixed_type2

    print("\n2. Contribution from reflections:")
    print(f"Sum = ({num_refl_type1}) * {num_colors}^({n//2}) + ({num_refl_type2}) * {num_colors}^({n//2}+1)")
    print(f"Sum for all {n} reflections = {reflection_sum}")
    
    total_sum_fixed_points = rotation_sum + reflection_sum
    print(f"\nTotal sum of fixed points = {rotation_sum} + {reflection_sum} = {total_sum_fixed_points}")

    num_orbits = total_sum_fixed_points // group_order
    
    print(f"\nNumber of non-isomorphic graphs = Total sum / {group_order}")
    print(f"Result = {total_sum_fixed_points} / {group_order} = {num_orbits}")


solve()
<<<538721474933832718133484343360293393333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333>>>