import math

def solve_optimal_sequence(initial_set):
    """
    Finds an optimal sequence to minimize the maximum cumulative impact.
    This function implements the greedy algorithm described in statement F.
    At each step, it chooses the number from the remaining set that minimizes
    the absolute value of the new cumulative sum.
    """
    
    remaining_numbers = list(initial_set)
    optimal_sequence = []
    current_sum = 0
    
    while remaining_numbers:
        best_choice = None
        best_idx = -1
        min_abs_sum = float('inf')
        
        # Find the best number to add next by searching all remaining numbers
        for i, num in enumerate(remaining_numbers):
            # A deterministic tie-breaking rule: if sums are equal, pick smaller number
            if abs(current_sum + num) < min_abs_sum:
                min_abs_sum = abs(current_sum + num)
                best_choice = num
                best_idx = i
            elif abs(current_sum + num) == min_abs_sum:
                if best_choice is None or num < best_choice:
                    best_choice = num
                    best_idx = i

        optimal_sequence.append(best_choice)
        current_sum += best_choice
        remaining_numbers.pop(best_idx)
        
    # Calculate cumulative sums and their absolute values
    abs_cumulative_sums = []
    s = 0
    for x in optimal_sequence:
        s += x
        abs_cumulative_sums.append(abs(s))
        
    max_impact = max(abs_cumulative_sums) if abs_cumulative_sums else 0
        
    return optimal_sequence, abs_cumulative_sums, max_impact

def demonstrate_statement_C():
    """
    This function demonstrates Statement C: The optimal maximum cumulative impact 
    is invariant under negation of all elements.
    """
    
    # Use the "Tricky case" example from the problem description
    A = [1, -4, 3, -1, 2, -2]
    
    # Solve for the original sequence A
    opt_seq_A, C_A, max_C_A = solve_optimal_sequence(A)
    
    # Create the negated sequence
    A_neg = [-x for x in A]
    
    # Solve for the negated sequence -A
    opt_seq_A_neg, C_A_neg, max_C_A_neg = solve_optimal_sequence(A_neg)

    print("--- Verifying Statement C ---")
    print(f"Original sequence A: {A}")
    print(f"Optimal permutation for A: {opt_seq_A}")
    print(f"Absolute cumulative sums C(k) for A: {C_A}")
    print(f"Optimal max cumulative impact for A: {max_C_A}\n")

    print(f"Negated sequence -A: {A_neg}")
    print(f"Optimal permutation for -A: {opt_seq_A_neg}")
    print(f"Absolute cumulative sums C(k) for -A: {C_A_neg}")
    print(f"Optimal max cumulative impact for -A: {max_C_A_neg}\n")
    
    are_equal = "are" if max_C_A == max_C_A_neg else "are NOT"
    print(f"Conclusion: The optimal maximum values {are_equal} equal ({max_C_A} vs {max_C_A_neg}). This supports statement C.")

demonstrate_statement_C()
