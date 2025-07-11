def solve_modulo_permutation(x, a):
    """
    Finds the permutation of a list 'a' that maximizes the final value of x
    after sequential modulo operations.

    The problem is to minimize |x_initial - x_final|, which is equivalent to
    maximizing x_final since x_final <= x_initial.

    This function uses dynamic programming. The state dp[mask] stores a dictionary
    mapping possible outcomes to the permutation of numbers (from the subset of 'a'
    represented by 'mask') that achieves that outcome.

    - dp: a list of dictionaries, where the index is the bitmask.
    - mask: an integer where the i-th bit is 1 if a[i] is in the subset.
    - dp[mask]: a dictionary {outcome: permutation_tuple}.

    The algorithm iterates through all subsets of 'a', from size 1 to n,
    calculating all possible outcomes for each subset.
    """
    n = len(a)
    if n == 0:
        print(f"Initial x: {x}")
        print("List 'a' is empty.")
        print(f"Final x: {x}")
        print("Equation: x = x")
        return

    dp = [{} for _ in range(1 << n)]

    # Base cases: masks with a single element
    for i in range(n):
        mask = 1 << i
        outcome = x % a[i]
        dp[mask][outcome] = (a[i],)

    # Iterate through masks from size 2 to n
    for i in range(1, 1 << n):
        for j in range(n):
            # If a[j] is in the current subset (mask i)
            if (i >> j) & 1:
                # Get the previous mask (subset without a[j])
                prev_mask = i ^ (1 << j)
                if not dp[prev_mask]:
                    continue
                
                # For each outcome from the previous subset,
                # calculate the new outcome by applying modulo a[j]
                for prev_outcome, prev_permutation in dp[prev_mask].items():
                    new_outcome = prev_outcome % a[j]
                    new_permutation = prev_permutation + (a[j],)
                    # Store the new outcome and its permutation
                    # If an outcome can be reached by multiple permutations,
                    # this will just keep the last one found.
                    dp[i][new_outcome] = new_permutation

    # The final mask represents the entire set 'a'
    full_mask = (1 << n) - 1
    final_results = dp[full_mask]

    if not final_results:
        best_x = x
        best_perm = []
    else:
        # Find the maximum possible final x
        best_x = max(final_results.keys())
        best_perm = final_results[best_x]

    # Output the results
    print(f"Initial x: {x}")
    print(f"List a: {a}")
    print(f"Best permutation found: {list(best_perm)}")
    print(f"Final x: {best_x}")

    equation_str = f"{x}"
    temp_x = x
    for val in best_perm:
        new_x = temp_x % val
        equation_str += f" % {val} (= {new_x})"
        temp_x = new_x
    
    print(f"Equation: {x} = {equation_str}")


# Example usage:
initial_x = 100
number_list = [10, 22, 23]
solve_modulo_permutation(initial_x, number_list)
