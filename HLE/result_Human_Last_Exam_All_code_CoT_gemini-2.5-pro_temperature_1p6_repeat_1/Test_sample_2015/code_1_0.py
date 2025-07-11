def find_cycle_and_fixed_values():
    """
    This function runs a simulation for all positive three-digit numbers to find the set of
    all values included in fixed points or cycles according to the specified process.

    The process is:
    1. For a number, create A (smallest digit permutation) and B (largest).
    2. The next number is B - A + 1.
    3. Repeat with the new number.
    """
    
    # This set will store all unique numbers found in any cycle or fixed point.
    final_cycle_values = set()

    # Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        
        current_num = start_num
        history = []
        
        # Loop until a number repeats in the sequence for this starting number.
        while current_num not in history:
            history.append(current_num)
            
            # Treat the number as a 3-digit string, padding with zeros if necessary.
            s_num = str(current_num).zfill(3)
            
            # Sort the digits to get the smallest and largest permutations.
            sorted_digits = sorted(list(s_num))
            
            # Calculate A (smallest value) and B (largest value).
            val_a = int("".join(sorted_digits))
            val_b = int("".join(reversed(sorted_digits)))
            
            # Calculate the next number in the sequence.
            current_num = val_b - val_a + 1
            
        # A cycle has been detected. The cycle begins at the first index of the repeated number.
        try:
            cycle_start_index = history.index(current_num)
            # The cycle consists of all elements from that index to the end of the history.
            cycle = history[cycle_start_index:]
            # Add all numbers from this cycle to our final set.
            final_cycle_values.update(cycle)
        except ValueError:
            # This block should not be reached due to the while loop's logic,
            # but is included for robustness.
            pass

    # Sort the final set of values in ascending order for the output.
    sorted_values = sorted(list(final_cycle_values))
    
    # Print the result in the specified format, e.g., {1, 2, 3}
    print(f"{{{', '.join(map(str, sorted_values))}}}")

find_cycle_and_fixed_values()