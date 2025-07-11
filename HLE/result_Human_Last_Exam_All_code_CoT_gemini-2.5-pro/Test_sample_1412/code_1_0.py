import math

def phi(n):
    """
    Computes Euler's Totient function phi(n), which counts the positive integers
    up to a given integer n that are relatively prime to n.
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
    Returns a list of all divisors of a number n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def solve():
    """
    Calculates the number of non-isomorphic graphs based on the plan.
    """
    n = 1000
    num_colors = 2

    # Part 1: Sum over rotations
    # The number of colorings fixed by a rotation of k positions in a cycle of length n
    # is num_colors^gcd(k, n).
    # Using Burnside's Lemma, the sum over all rotations is sum_{k=0 to n-1} num_colors^gcd(k, n).
    # This can be more efficiently calculated as sum_{d|n} phi(n/d) * num_colors^d.
    divisors = get_divisors(n)
    sum_rotations = 0
    for d in divisors:
        term = phi(n // d) * (num_colors**d)
        sum_rotations += term

    # Part 2: Sum over reflections
    # For an even n, there are two types of reflections in the dihedral group D_n.
    # n/2 reflections through opposite vertices: n/2 * num_colors^(n/2)
    # n/2 reflections through midpoints of opposite edges: n/2 * num_colors^(n/2 + 1)
    sum_reflections_vertices = (n // 2) * (num_colors**(n // 2))
    sum_reflections_edges = (n // 2) * (num_colors**(n // 2 + 1))
    sum_reflections = sum_reflections_vertices + sum_reflections_edges

    # Part 3: Apply Burnside's Lemma
    # Total sum of fixed points
    total_sum_fixed_points = sum_rotations + sum_reflections
    
    # The size of the dihedral group D_n is 2n.
    group_size = 2 * n
    
    # The number of non-isomorphic configurations is the average number of fixed points.
    num_non_isomorphic = total_sum_fixed_points // group_size

    # Print the components of the final calculation as requested
    print("This problem is equivalent to counting non-isomorphic 2-edge-colorings of a cycle with 1000 vertices.")
    print("Using Burnside's Lemma for the Dihedral Group D_1000:")
    print(f"Number of vertices (n): {n}")
    print(f"Order of the group D_n (|G|): {group_size}")
    print(f"Sum of fixed points for rotations: {sum_rotations}")
    print(f"Sum of fixed points for reflections: {sum_reflections}")
    print("\nFinal Calculation:")
    print(f"({sum_rotations} + {sum_reflections}) / {group_size} = {num_non_isomorphic}")
    
    # The final answer in the required format
    print("\nFinal Answer:")
    print(f"<<<{num_non_isomorphic}>>>")

solve()