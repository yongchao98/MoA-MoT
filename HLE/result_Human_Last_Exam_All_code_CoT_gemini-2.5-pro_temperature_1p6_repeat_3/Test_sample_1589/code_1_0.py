import itertools

def solve(n):
    """
    Calculates a(n), the number of valid pairings.

    A pairing is a permutation y of {n+1, ..., 2n} such that for the
    pairs (i, y_i), the 2n values {y_i + i} and {y_i - i} are all distinct.
    """
    count = 0
    
    # The set of numbers to be paired with {1, ..., n}
    y_set = range(n + 1, 2 * n + 1)
    
    # Iterate through all possible permutations of y_set
    # Each permutation 'p' corresponds to a possible pairing (1,p[0]), (2,p[1]),...
    for p in itertools.permutations(y_set):
        
        # A set to store all the generated values to check for uniqueness
        all_values = set()
        is_valid = True
        
        # For each number from 1 to n, generate the sum and difference
        for i in range(n):
            # x_i is the number from the first set, {1, ..., n}
            x_i = i + 1
            # y_i is the corresponding number from the permutation of {n+1, ..., 2n}
            y_i = p[i]
            
            s = y_i + x_i  # The sum
            d = y_i - x_i  # The difference
            
            # Check if either the sum or difference has been seen before
            if s in all_values or d in all_values:
                is_valid = False
                break
            
            # If they are new, add them to the set of seen values
            all_values.add(s)
            all_values.add(d)
        
        # If the inner loop completed without finding duplicates, this permutation is valid
        if is_valid:
            count += 1
            
    return count

if __name__ == '__main__':
    n = 10
    result = solve(n)
    print(f"a({n}) = {result}")
