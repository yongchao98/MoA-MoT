def solve_digit_cycle_problem():
    """
    This function simulates the described digit reordering process for all
    three-digit numbers to find the set of all values that are either
    fixed points or part of a cycle.
    """
    final_cycle_values = set()

    # Iterate through every positive three-digit number.
    for start_num in range(100, 1000):
        current_num = start_num
        history = []

        # Keep generating new numbers until we find one that's already in the
        # history for this specific starting number.
        while current_num not in history:
            history.append(current_num)

            # To handle numbers with fewer than 3 digits (like 1, 99),
            # we format them as a 3-digit string with leading zeros.
            s_num = str(current_num).zfill(3)
            
            # Sort the digits to get the smallest number A.
            sorted_digits = sorted(list(s_num))
            A = int("".join(sorted_digits))
            
            # Reverse the sorted digits to get the largest number B.
            B = int("".join(reversed(sorted_digits)))
            
            # Calculate the next number in the sequence.
            current_num = B - A + 1

        # A cycle or fixed point is found. The loop has terminated because
        # current_num is already in history.
        # We find the first time it appeared to identify the start of the cycle.
        cycle_start_index = history.index(current_num)
        
        # The cycle consists of all numbers from that first appearance to the end.
        cycle = history[cycle_start_index:]
        
        # Add all numbers from this cycle to our final set.
        # The set automatically handles duplicates.
        final_cycle_values.update(cycle)

    # Sort the final set of numbers in ascending order for the output.
    sorted_values = sorted(list(final_cycle_values))
    
    # Format the output string as requested, e.g., {1, 2, 3}
    result_string = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_string)

solve_digit_cycle_problem()