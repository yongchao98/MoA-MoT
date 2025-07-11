def find_cycle_values():
    """
    Finds all values that are part of a fixed point or cycle for the B-A+1 process 
    on all positive 3-digit numbers.
    """
    
    # This set will store all numbers that are part of any cycle or fixed point.
    cycle_values = set()

    # Iterate through all positive 3-digit numbers.
    for i in range(100, 1000):
        
        # If the starting number leads to a path we've already explored, we can skip it.
        # This check is not strictly necessary but is an optimization. We can check within the loop.

        path_history = []
        current_num = i
        
        # Generate the sequence for the current starting number.
        while True:
            # If the current number is in our master set, its cycle is already known.
            # We can stop exploring this path.
            if current_num in cycle_values:
                break
            
            # If we've seen the current number before in this specific path, we've found a new cycle.
            if current_num in path_history:
                # Identify the start of the cycle.
                cycle_start_index = path_history.index(current_num)
                # Extract the members of the cycle.
                new_cycle_members = path_history[cycle_start_index:]
                # Add these members to our master set.
                cycle_values.update(new_cycle_members)
                break

            # Add the current number to this path's history.
            path_history.append(current_num)

            # Pad the number to 3 digits to handle cases like '1' -> '001'.
            s_num = str(current_num).zfill(3)
            
            # Sort digits to get A and B.
            sorted_digits = sorted(list(s_num))
            
            val_a = int("".join(sorted_digits))
            val_b = int("".join(reversed(sorted_digits)))

            # Calculate the next number in the sequence.
            current_num = val_b - val_a + 1

    # Sort the final set of values in ascending order.
    sorted_result = sorted(list(cycle_values))
    
    # Format the output as requested, e.g., {1, 2, 3}.
    print(f"{{{', '.join(map(str, sorted_result))}}}")

if __name__ == '__main__':
    find_cycle_values()