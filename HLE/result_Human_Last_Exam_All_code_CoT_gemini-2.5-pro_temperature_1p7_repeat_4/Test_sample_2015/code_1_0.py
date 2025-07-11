def find_cycle_values():
    """
    This function calculates the total set of values included in fixed points or cycles
    for a specific process run on all positive three-digit numbers.
    """
    # A set to store all numbers that are found to be part of a cycle or are fixed points.
    cycle_and_fixed_values = set()

    # Iterate through all possible three-digit starting numbers.
    for start_num in range(100, 1000):
        
        # Optimization: If we've already processed this number as part of a previous path, skip it.
        if start_num in cycle_and_fixed_values:
            continue
            
        current_num = start_num
        # This list tracks the path taken from the current starting number.
        path = []

        # Continue the process until a number is repeated or we hit a known cycle.
        while current_num not in path and current_num not in cycle_and_fixed_values:
            path.append(current_num)

            # Ensure the number is treated as having 3 digits by padding with leading zeros.
            s_num = str(current_num).zfill(3)
            
            # Sort the digits to get the characters for A and B.
            digits = sorted(list(s_num))
            
            # Create A, the smallest number, and B, the largest number.
            a = int("".join(digits))
            b = int("".join(reversed(digits)))
            
            # Calculate the next number in the sequence.
            current_num = b - a + 1
        
        # If the loop terminated because a number was found in the current path,
        # it means we have discovered a new cycle.
        if current_num in path:
            # Find the starting index of the cycle.
            cycle_start_index = path.index(current_num)
            # Extract the numbers that form the cycle.
            newly_found_cycle_members = path[cycle_start_index:]
            # Add these newly found numbers to our master set.
            cycle_and_fixed_values.update(newly_found_cycle_members)

    # Convert the set to a list and sort it in ascending order for the final output.
    sorted_values = sorted(list(cycle_and_fixed_values))
    
    # Format and print the result as a set literal.
    print(f"{{{', '.join(map(str, sorted_values))}}}")

find_cycle_values()