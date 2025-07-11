import itertools

def find_optimal_sequence(initial_list):
    """
    Finds the optimal permutation of a list of numbers to minimize the
    maximum absolute cumulative sum. This implementation uses brute force
    and is suitable only for small lists.

    Args:
        initial_list: A list of numbers.

    Returns:
        A tuple containing the optimal permutation, the sequence of
        absolute cumulative sums, and the minimized maximum value.
    """
    if not initial_list:
        return [], [], 0

    min_max_abs_sum = float('inf')
    optimal_permutation = None
    optimal_abs_cum_sums = None

    # Iterate through all unique permutations of the input list
    for p in set(itertools.permutations(initial_list)):
        cum_sum = 0
        max_abs_sum = 0
        current_abs_cum_sums = []
        
        for x in p:
            cum_sum += x
            abs_val = abs(cum_sum)
            current_abs_cum_sums.append(abs_val)
            if abs_val > max_abs_sum:
                max_abs_sum = abs_val

        if max_abs_sum < min_max_abs_sum:
            min_max_abs_sum = max_abs_sum
            optimal_permutation = p
            optimal_abs_cum_sums = current_abs_cum_sums

    return optimal_permutation, optimal_abs_cum_sums, min_max_abs_sum

# --- Example Execution ---
# We use the "Tricky case" from the prompt to demonstrate the code.
# Note: The examples in the prompt contain some calculation errors.
# This code will find the correct optimal solution.
input_sequence = [1, -4, 3, -1, 2, -2]

opt_perm, opt_sums, min_max = find_optimal_sequence(input_sequence)

print(f"Input Sequence: {input_sequence}")
print("-" * 30)
print(f"Optimal Permutation: {list(opt_perm)}")
# The following line demonstrates the "output each number" instruction by showing
# the sequence of absolute cumulative sums that results from the optimal permutation.
print(f"Absolute Cumulative Sums: {opt_sums}")
print(f"Minimized Maximum Value: {min_max}")
