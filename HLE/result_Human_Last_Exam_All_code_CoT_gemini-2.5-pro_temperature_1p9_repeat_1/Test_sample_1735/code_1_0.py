import math

def solve_optimal_sequence(numbers):
    """
    Finds a near-optimal permutation of a sequence to minimize the maximum
    absolute cumulative sum using a greedy algorithm.
    At each step, it selects the number that minimizes the absolute value of the
    new cumulative sum.

    Args:
        numbers (list): A list of real numbers.

    Returns:
        None. Prints the resulting permutation and its cumulative sums.
    """
    
    # Copy the input list to modify it
    remaining_numbers = list(numbers)
    
    permutation = []
    cumulative_sum = 0
    cumulative_sums_abs = []

    for _ in range(len(numbers)):
        best_num = None
        min_abs_sum = float('inf')

        # Find the number that makes the next cumulative sum closest to zero
        for num in remaining_numbers:
            if abs(cumulative_sum + num) < min_abs_sum:
                min_abs_sum = abs(cumulative_sum + num)
                best_num = num
            # Tie-breaking rule: prefer smaller magnitude numbers
            elif abs(cumulative_sum + num) == min_abs_sum:
                if abs(num) < abs(best_num):
                    best_num = num

        # Update the state
        cumulative_sum += best_num
        permutation.append(best_num)
        cumulative_sums_abs.append(abs(cumulative_sum))
        remaining_numbers.remove(best_num)

    print(f"Input: {numbers}")
    print(f"Optimal: {permutation}")
    print(f"Absolute Cumulative Sums: {cumulative_sums_abs}")
    print(f"Maximum absolute cumulative sum: {max(cumulative_sums_abs)}")


# Edge case from the problem description
input_sequence = [2, 2, -3, -3, 2]
solve_optimal_sequence(input_sequence)