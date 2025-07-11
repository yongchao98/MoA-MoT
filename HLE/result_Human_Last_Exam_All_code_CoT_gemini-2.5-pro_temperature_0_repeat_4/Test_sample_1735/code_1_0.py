import math

def solve_dynamic_sequence(numbers):
    """
    Finds an optimal permutation of a sequence of numbers to minimize the
    maximum absolute cumulative sum, using a greedy algorithm.

    Args:
        numbers: A list of real numbers.
    """
    # Use a copy of the list to modify it during the process
    remaining_numbers = list(numbers)
    
    solution_sequence = []
    cumulative_sums = []
    abs_cumulative_sums = []
    current_sum = 0

    print(f"Input: {numbers}")

    while remaining_numbers:
        best_next_num = None
        min_abs_sum = float('inf')

        # Find the number that minimizes the next absolute cumulative sum
        for num in remaining_numbers:
            abs_sum = abs(current_sum + num)
            if abs_sum < min_abs_sum:
                min_abs_sum = abs_sum
                best_next_num = num
            # Simple tie-breaking: if sums are equal, prefer smaller numbers
            elif abs_sum == min_abs_sum:
                if best_next_num is None or num < best_next_num:
                    best_next_num = num

        # Add the best number to the solution and update state
        solution_sequence.append(best_next_num)
        current_sum += best_next_num
        cumulative_sums.append(current_sum)
        abs_cumulative_sums.append(abs(current_sum))
        remaining_numbers.remove(best_next_num)

    # Output the results as per the example format
    # The problem asks to output each number in the final equation.
    # We will print the optimal sequence and the resulting cumulative values.
    
    # Format the optimal sequence string like {-1, 2, -4, 3}
    optimal_str = "{" + ", ".join(map(str, solution_sequence)) + "}"
    print(f"Optimal reordering: {optimal_str}")

    # Format the absolute cumulative sums string like {1, 1, 3, 0}
    abs_sums_str = "{" + ", ".join(map(str, abs_cumulative_sums)) + "}"
    print(f"Absolute cumulative sums: {abs_sums_str}")
    
    print(f"Minimized maximum value: {max(abs_cumulative_sums)}")


# Example from the problem description
input_sequence = [3, -1, -4, 2]
solve_dynamic_sequence(input_sequence)
