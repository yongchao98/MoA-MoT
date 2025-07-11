import math

def solve_optimal_sequence(numbers):
    """
    Finds an optimal permutation of a sequence of numbers to minimize the
    maximum absolute cumulative sum, using a greedy algorithm.
    """
    
    remaining_numbers = list(numbers)
    optimal_sequence = []
    cumulative_sums = []
    abs_cumulative_sums = []
    current_sum = 0

    print(f"Input: {numbers}")

    for _ in range(len(numbers)):
        best_choice = None
        min_abs_sum = float('inf')

        # Find the number in the remaining list that minimizes the next cumulative sum
        for num in remaining_numbers:
            if abs(current_sum + num) < min_abs_sum:
                min_abs_sum = abs(current_sum + num)
                best_choice = num
            # Tie-breaking: if sums are equal, pick the smaller number
            elif abs(current_sum + num) == min_abs_sum:
                if best_choice is None or num < best_choice:
                    best_choice = num

        # Add the chosen number to our sequence
        optimal_sequence.append(best_choice)
        current_sum += best_choice
        cumulative_sums.append(current_sum)
        abs_cumulative_sums.append(abs(current_sum))
        
        # Remove it from the list of available numbers
        remaining_numbers.remove(best_choice)

    print(f"Optimal reordering: {optimal_sequence}")
    print(f"Cumulative sums: {cumulative_sums}")
    print(f"Absolute cumulative sums: {abs_cumulative_sums}")
    max_impact = max(abs_cumulative_sums) if abs_cumulative_sums else 0
    print(f"The minimum maximum cumulative impact is: {max_impact}")

# Example from the problem description: Tricky case
tricky_case = [1, -4, 3, -1, 2, -2]
solve_optimal_sequence(tricky_case)

<<<A, C, F, H>>>