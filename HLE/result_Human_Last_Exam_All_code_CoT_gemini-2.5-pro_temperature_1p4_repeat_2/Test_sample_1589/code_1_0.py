import itertools

def solve_a10():
    """
    This function calculates a(10), which is the number of ways to pair 
    the natural numbers from 1 to 10 with those from 11 to 20 into 
    n pairs (i, y_i) such that the 2n numbers y_i+i and y_i-i are all different.
    
    The function works by checking every possible permutation.
    Note: This computation may take a few moments as it needs to check
    10! = 3,628,800 permutations.
    """
    n = 10
    count = 0
    
    # The set of numbers y_i are from {n+1, ..., 2n}
    y_values = range(n + 1, 2 * n + 1)
    
    # The set of indices i are from {1, ..., n}
    x_values = range(1, n + 1)
    
    # Generate an iterator for all permutations of y_values.
    all_permutations = itertools.permutations(y_values)
    
    # Iterate through each permutation p = (y_1, y_2, ..., y_n)
    for p in all_permutations:
        generated_values = set()
        
        # For each index i, we have a pair (i, y_i).
        # We calculate the two required values and add them to a set.
        for i in range(n):
            x = x_values[i]  # This corresponds to i=1, 2, ..., 10
            y = p[i]         # This is the corresponding y_i
            
            generated_values.add(y + x)
            generated_values.add(y - x)
            
        # The condition is that all 2n generated values must be different.
        # This is true if the size of the set is equal to 2n.
        if len(generated_values) == 2 * n:
            count += 1
            
    # Print the final result in the format: a(10) = result
    print(f"a(10) = {count}")

# Run the solver function
solve_a10()