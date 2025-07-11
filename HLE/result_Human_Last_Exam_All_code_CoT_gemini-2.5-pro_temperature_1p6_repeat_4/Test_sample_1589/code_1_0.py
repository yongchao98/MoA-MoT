import itertools

def solve_pairing_problem():
    """
    Calculates a(n) for n=10.
    a(n) is the number of ways to pair the natural numbers from 1 to n
    with those between n+1 and 2n into n pairs (i, y_i) such that
    the 2n numbers y_i+i and y_i-i are all different.
    """
    n = 10
    count = 0
    
    # The set of numbers {y_i} is a permutation of {n+1, ..., 2n}
    y_values_to_permute = range(n + 1, 2 * n + 1)
    
    # The set of numbers {i} is {1, ..., n}
    i_values = range(1, n + 1)
    
    # Create an iterator for all permutations of y_values
    all_permutations = itertools.permutations(y_values_to_permute)
    
    # Loop through each permutation
    # p will be a tuple like (y_1, y_2, ..., y_10)
    for p in all_permutations:
        
        generated_values = set()
        is_valid_permutation = True
        
        # For each pair (i, y_i), calculate the two required values
        for i, y_i in zip(i_values, p):
            val_plus = y_i + i
            val_minus = y_i - i
            
            # Check for collisions. If a value is already in the set,
            # this permutation is invalid.
            if val_plus in generated_values:
                is_valid_permutation = False
                break
            generated_values.add(val_plus)
            
            if val_minus in generated_values:
                is_valid_permutation = False
                break
            generated_values.add(val_minus)
        
        # If the inner loop completed without finding a collision,
        # then this permutation is valid.
        if is_valid_permutation:
            count += 1
            
    # As requested, output the final equation with the numbers.
    print(f"a({n}) = {count}")

# Execute the function to find and print the solution.
# Note: This script will take a few moments to run due to the large
# number of permutations to check (10! = 3,628,800).
solve_pairing_problem()