import collections

def find_cycle_values():
    """
    This function executes the described process for all three-digit numbers
    to find the complete set of values that are part of cycles or are fixed points.
    """
    
    # A set to store all numbers that are part of a cycle or are fixed points.
    cycle_and_fixed_values = set()

    # An optimization to keep track of numbers whose paths have already been fully explored.
    visited_numbers = set()

    # Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        # If we have already processed this number as part of a previous path, skip it.
        if start_num in visited_numbers:
            continue
        
        # Store the sequence of numbers for the current starting number.
        path = []
        current_num = start_num
        
        # Trace the sequence until we find a number we've seen before in this specific path.
        while current_num not in path:
            # Optimization: if we hit a number that is part of a known cycle,
            # we can stop this path, as we know it will lead to that same cycle.
            if current_num in cycle_and_fixed_values:
                break
                
            path.append(current_num)
            
            # Calculate the next number in the sequence.
            # We treat all numbers as 3-digit numbers by padding with zeros.
            s_num = str(current_num).zfill(3)
            
            # Get the digits and sort them.
            digits = sorted(list(s_num))
            
            # Form A (smallest value) and B (largest value).
            A = int("".join(digits))
            B = int("".join(reversed(digits)))
            
            # Calculate the next number in the sequence.
            current_num = B - A + 1
        else:
            # This block executes if the while loop completes without a 'break'.
            # This means a new cycle was found because 'current_num' is already in 'path'.
            try:
                # Find where the cycle begins in our path.
                cycle_start_index = path.index(current_num)
                # Extract all the numbers that form the cycle.
                new_cycle_members = path[cycle_start_index:]
                # Add these newly found cycle members to our global set.
                cycle_and_fixed_values.update(new_cycle_members)
            except ValueError:
                # This case should not be reached with the current logic.
                pass

        # Mark all numbers in the path we just traced as visited to avoid re-computation.
        visited_numbers.update(path)

    # Sort the final set of values in ascending order.
    final_values = sorted(list(cycle_and_fixed_values))
    
    # Format the output string as requested.
    # This prints each number in the final set.
    result_string = "{" + ", ".join(map(str, final_values)) + "}"
    print(result_string)

find_cycle_values()