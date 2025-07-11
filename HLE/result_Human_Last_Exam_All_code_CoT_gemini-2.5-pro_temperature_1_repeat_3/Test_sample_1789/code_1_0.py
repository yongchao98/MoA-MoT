def solve_antichain_union():
    """
    Calculates the size of the largest union of 20 antichains in
    [1, 823564528378596] in the divisor poset.

    This is equivalent to finding the number of integers m in the range
    such that the exponent of 2 in the prime factorization of m is less than 20.
    """
    n = 823564528378596
    k = 20

    total_size = 0
    terms = []
    
    # We sum the number of odd integers up to n / 2^j for j from 0 to k-1.
    for j in range(k):
        power_of_2 = 1 << j  # Efficient way to compute 2**j
        
        limit = n // power_of_2
        
        # The number of odd integers in [1, limit] is (limit + 1) // 2
        num_odds = (limit + 1) // 2
        
        terms.append(num_odds)
        total_size += num_odds
        
    # Format the equation string with each term
    equation_str = " + ".join(map(str, terms))
    
    print("The size is calculated by the following sum:")
    print(f"{equation_str} = {total_size}")

solve_antichain_union()