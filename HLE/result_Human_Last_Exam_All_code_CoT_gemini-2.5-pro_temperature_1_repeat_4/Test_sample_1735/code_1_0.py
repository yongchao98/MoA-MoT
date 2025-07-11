import math

def calculate_max_abs_cumulative_sum(seq):
    """Calculates and prints the max absolute cumulative sum for a sequence."""
    if not seq:
        return 0
    current_sum = 0
    max_abs_sum = 0
    cumulative_sums = []
    print(f"Sequence: {seq}")
    
    # Build the string showing the calculation for clarity
    calc_str = "Cumulative sums: "
    for x in seq:
        current_sum += x
        cumulative_sums.append(current_sum)
        calc_str += f"{current_sum} "
        if abs(current_sum) > max_abs_sum:
            max_abs_sum = abs(current_sum)
            
    print(calc_str)
    print(f"Maximum Absolute Cumulative Sum: {max_abs_sum}\n")
    return max_abs_sum

def simple_greedy_solver(numbers):
    """
    Implements the simple greedy strategy from statement F:
    At each step k, choose the next number to minimize |current_sum + next_number|.
    """
    remaining = list(numbers)
    solution = []
    current_sum = 0
    
    while remaining:
        best_choice = None
        # Initialize min_abs_sum with a very large number
        min_abs_sum = float('inf')
        
        for num in remaining:
            abs_sum = abs(current_sum + num)
            if abs_sum < min_abs_sum:
                min_abs_sum = abs_sum
                best_choice = num
            # A simple tie-breaking rule (e.g., choose smaller number) could be added
            # but is not necessary for this counterexample.
        
        solution.append(best_choice)
        remaining.remove(best_choice)
        current_sum += best_choice
        
    return solution

# --- Analysis of Statement F ---
print("--- Analysis of Statement F ---")
print("Statement F claims a simple greedy algorithm is always optimal.")
print("We test this with the counterexample set {10, 8, -15, -9}.\n")

input_set = {10, 8, -15, -9}

# Run the simple greedy algorithm described in statement F
print("1. Running the simple greedy algorithm (from statement F):")
greedy_solution = simple_greedy_solver(input_set)
greedy_max_sum = calculate_max_abs_cumulative_sum(greedy_solution)

# Test the known better solution
print("2. Testing a known better permutation:")
better_solution = [8, -15, 10, -9]
better_max_sum = calculate_max_abs_cumulative_sum(better_solution)

# Conclusion
if better_max_sum < greedy_max_sum:
    print("--- Conclusion ---")
    print(f"The simple greedy algorithm found a solution with max sum {greedy_max_sum},")
    print(f"but a different permutation exists with a smaller max sum of {better_max_sum}.")
    print("This demonstrates that the greedy strategy of minimizing C(k) at each step is not always optimal.")
    print("Therefore, statement F is FALSE.")
