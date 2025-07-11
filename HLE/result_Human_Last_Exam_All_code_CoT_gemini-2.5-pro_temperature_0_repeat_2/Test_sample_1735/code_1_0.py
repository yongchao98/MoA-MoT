import math

def solve_optimizing_dynamic_sequences(A):
    """
    Finds an optimal permutation of a sequence A to minimize the maximum
    absolute cumulative sum, using a greedy algorithm.
    """
    # Separate into positive and negative numbers
    # Sort them by descending absolute value
    pos = sorted([x for x in A if x >= 0], reverse=True)
    neg = sorted([x for x in A if x < 0])

    optimal_sequence = []
    current_sum = 0

    # Build the optimal sequence greedily
    while pos or neg:
        # If current_sum is non-negative, we prefer to add a negative number.
        # If current_sum is negative, we prefer to add a positive number.
        use_neg = (current_sum >= 0)

        # Handle edge cases where one list is empty
        if not neg:
            use_neg = False
        if not pos:
            use_neg = True
            
        # Tie-breaking for current_sum == 0
        if current_sum == 0 and pos and neg:
            if pos[0] < abs(neg[0]):
                use_neg = False
            else:
                use_neg = True

        if use_neg:
            num = neg.pop(0)
        else:
            num = pos.pop(0)
        
        optimal_sequence.append(num)
        current_sum += num

    # Print the results
    print(f"Input: {A}")
    print(f"Optimal sequence: {optimal_sequence}")
    print("Cumulative sums:")
    
    s = 0
    max_abs_sum = 0
    sum_str_list = []
    for x in optimal_sequence:
        s += x
        if x < 0:
            sum_str_list.append(f"- {abs(x)}")
        else:
            sum_str_list.append(f"+ {x}")
        
        # Join the elements for the equation string, handling the first element
        if len(sum_str_list) == 1:
            equation_str = f"{optimal_sequence[0]}"
        else:
            # Start with the first number, then add operators for subsequent numbers
            equation_str = f"{optimal_sequence[0]} "
            for part in sum_str_list[1:]:
                equation_str += f"{part} "
        
        print(f"| {equation_str}| = {abs(s)}")
        if abs(s) > max_abs_sum:
            max_abs_sum = abs(s)
            
    print(f"Maximum absolute cumulative sum: {max_abs_sum}")

# Tricky case from the problem description
tricky_case_input = [1, -4, 3, -1, 2, -2]
solve_optimizing_dynamic_sequences(tricky_case_input)