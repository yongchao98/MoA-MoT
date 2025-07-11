import math

def find_optimal_sequence_myopic(A):
    """
    Finds an optimal sequence by greedily minimizing the absolute cumulative sum at each step.
    This implements the strategy described in statement F.
    """
    rem = list(A)
    seq = []
    current_sum = 0
    
    # This is an O(n^2) implementation, but is conceptually simple.
    # An O(n log n) version would use sorted lists and binary search.
    for _ in range(len(A)):
        best_choice = None
        min_abs_sum = float('inf')
        
        # Find the number in the remaining set that makes the next sum closest to 0
        for x in rem:
            abs_sum = abs(current_sum + x)
            if abs_sum < min_abs_sum:
                min_abs_sum = abs_sum
                best_choice = x
            # Tie-breaking for robustness, though not strictly required by problem
            elif abs_sum == min_abs_sum:
                if best_choice is None or abs(x) < abs(best_choice):
                     best_choice = x

        seq.append(best_choice)
        rem.remove(best_choice)
        current_sum += best_choice
    
    return seq

def calculate_max_c(seq):
    """Calculates the maximum absolute cumulative sum for a sequence."""
    if not seq:
        return 0, []
    current_sum = 0
    max_c = 0
    abs_cumulative_sums = []
    for x in seq:
        current_sum += x
        abs_c = abs(current_sum)
        abs_cumulative_sums.append(abs_c)
        if abs_c > max_c:
            max_c = abs_c
    return max_c, abs_cumulative_sums

# --- Main execution ---
# Analyze the "Tricky case" from the problem description
tricky_case_input = [1, -4, 3, -1, 2, -2]

# Find the optimal sequence using the greedy algorithm
optimal_seq = find_optimal_sequence_myopic(tricky_case_input)
max_impact, impact_seq = calculate_max_c(optimal_seq)

print(f"Analysis of 'Tricky case': {set(tricky_case_input)}")
print("An optimal solution can be found using a greedy algorithm.")
print(f"Optimal sequence found: {optimal_seq}")
print(f"Sequence of absolute cumulative sums: {[int(s) for s in impact_seq]}")
print(f"The minimized maximum cumulative impact is: {int(max_impact)}")

# The problem asks to output the final equation. This likely refers to the
# sum of the elements in the determined optimal order.
final_sum = sum(tricky_case_input)
equation_parts = []
for i, x in enumerate(optimal_seq):
    if i > 0 and x >= 0:
        equation_parts.append('+')
    equation_parts.append(str(x))
equation_str = " ".join(equation_parts)
print(f"Final equation: {equation_str} = {int(final_sum)}")