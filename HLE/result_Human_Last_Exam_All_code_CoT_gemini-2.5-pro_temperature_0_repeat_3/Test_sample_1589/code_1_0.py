import itertools

def solve_a10():
    """
    This function calculates a(n) for n=10.
    a(n) is the number of ways to pair the numbers {1, ..., n} with {n+1, ..., 2n}
    such that for each pair (i, y_i), the 2n values {y_i + i, y_i - i} are all distinct.

    This problem is equivalent to finding the number of permutations p of {0, ..., n-1}
    such that the values {p[i] + i} are all distinct and {p[i] - i} are all distinct.
    """
    n = 10
    count = 0
    
    # We will iterate through all permutations of (0, 1, ..., n-1).
    # Using 0-based indexing is more convenient in Python.
    nums = range(n)
    
    # Generate all n! permutations.
    for p in itertools.permutations(nums):
        # p is a tuple representing a permutation, e.g., (1, 0, 2, ...)
        
        # Condition 1: The values p[i] - i must all be distinct.
        # We use a set to check for uniqueness efficiently.
        diffs = {p[i] - i for i in range(n)}
        
        # If the number of unique differences is less than n, this permutation is invalid.
        if len(diffs) < n:
            continue
            
        # Condition 2: The values p[i] + i must all be distinct.
        sums = {p[i] + i for i in range(n)}
        
        # If the number of unique sums is less than n, this permutation is invalid.
        if len(sums) < n:
            continue
            
        # If both conditions are met, we have found a valid pairing.
        count += 1
        
    # The problem asks to output the final equation.
    print(f"a({n}) = {count}")

# Execute the function to find and print the solution.
solve_a10()