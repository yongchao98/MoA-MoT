import itertools

def solve_and_demonstrate(numbers):
    """
    Solves the optimizing dynamic sequences problem using brute-force and demonstrates the result.
    This approach is only feasible for small n due to its O(n!) complexity.
    """
    if not numbers:
        print("Input sequence is empty.")
        return

    min_max_c = float('inf')
    optimal_sequence = None
    
    # Use set to handle duplicate permutations automatically
    unique_permutations = set(itertools.permutations(numbers))

    # Iterate through all unique permutations to find the optimal one
    for p in unique_permutations:
        max_c = 0
        current_sum = 0
        for x in p:
            current_sum += x
            if abs(current_sum) > max_c:
                max_c = abs(current_sum)
        
        if max_c < min_max_c:
            min_max_c = max_c
            optimal_sequence = p

    print(f"Input sequence: {numbers}")
    print(f"Optimal sequence found: {list(optimal_sequence)}")
    print(f"Minimum maximum cumulative impact: {min_max_c}")
    print("-" * 20)
    print("Calculation of cumulative sums for the optimal sequence:")

    # Print the step-by-step "equation" for the cumulative sum
    current_sum = 0
    for i, x in enumerate(optimal_sequence):
        prev_sum = current_sum
        current_sum += x
        # Use a consistent sign for the operation
        op = "+" if x >= 0 else "-"
        abs_x = abs(x)
        print(f"Step {i+1}: {prev_sum} {op} {abs_x} = {current_sum}")

# Using the first example from the problem description
# The prompt's sequence is {3, -1, -4, 2}
# The prompt's suggested "optimal" sequence is {-1, 2, -4, 3} which yields a max cumulative impact of 3.
# This script will find the true optimal value.
input_numbers = [3, -1, -4, 2]
solve_and_demonstrate(input_numbers)