def find_cycle_values():
    """
    This function runs the described process for all three-digit numbers to find
    the total set of values included in fixed points or cycles.
    """
    # This set will store all unique numbers that are part of a cycle or are a fixed point.
    all_cycle_values = set()

    # Iterate through every possible three-digit starting number.
    for start_num in range(100, 1000):
        
        current_num = start_num
        path_history = []
        
        # Keep track of the sequence to detect cycles.
        # We use a list to easily extract the cycle members once detected.
        while current_num not in path_history:
            path_history.append(current_num)
            
            # Treat the number as a 3-digit number by padding with leading zeros if necessary.
            # This is crucial for correctly sorting the digits. e.g., 99 -> "099".
            s_num = str(current_num).zfill(3)
            
            # Get the digits and sort them to form A and B.
            digits = sorted(list(s_num))
            
            # A is the number with digits in ascending order.
            A = int("".join(digits))
            
            # B is the number with digits in descending order.
            B = int("".join(digits[::-1]))
            
            # Calculate the next value in the sequence.
            current_num = B - A + 1
            
        # Once a number repeats, a cycle is found. 'current_num' is the first repeated value.
        # Find where the cycle begins in our path history.
        cycle_start_index = path_history.index(current_num)
        
        # The members of the cycle are all numbers from the start of the cycle to the end of the path.
        members_of_cycle = path_history[cycle_start_index:]
        
        # Add all members of the found cycle/fixed point to our master set.
        all_cycle_values.update(members_of_cycle)

    # Sort the final set of values in ascending order.
    sorted_values = sorted(list(all_cycle_values))
    
    # Format the output as a string representing a set, e.g., {1, 2, 3}
    result_string = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_string)

find_cycle_values()