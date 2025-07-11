import collections

def find_cycle_values():
    """
    This function simulates the described number process for all 3-digit numbers
    to find the set of all values that are part of a cycle or are a fixed point.
    """
    
    # This set will store the final collection of all numbers that are part of a cycle.
    final_cycle_members = set()

    def get_next_value(n):
        """
        Calculates the next value in the sequence for a given number n.
        The number is treated as a 3-digit number.
        """
        # Format the number as a 3-digit string, padding with leading zeros if needed.
        s_n = str(n).zfill(3)
        
        # Sort the digits to get the smallest number A
        sorted_digits = sorted(list(s_n))
        s_a = "".join(sorted_digits)
        a = int(s_a)
        
        # Reverse the sorted digits to get the largest number B
        s_b = "".join(sorted_digits[::-1])
        b = int(s_b)
        
        # Calculate and return the new value
        return b - a + 1

    # Iterate through every positive three-digit number.
    for i in range(100, 1000):
        path = []
        current_num = i
        
        # Follow the path until we encounter a number we've already seen in this specific path.
        while current_num not in path:
            path.append(current_num)
            current_num = get_next_value(current_num)
            
        # Once a number repeats, a cycle has been found.
        # We need to identify the numbers that form the cycle.
        try:
            cycle_start_index = path.index(current_num)
            cycle = path[cycle_start_index:]
            
            # Add all numbers in this newly found cycle to our master set.
            final_cycle_members.update(cycle)
        except ValueError:
            # This case should not be reached with the current logic.
            pass

    # Sort the collected numbers in ascending order.
    sorted_members = sorted(list(final_cycle_members))
    
    # Format the output string as requested.
    result_str = "{" + ", ".join(map(str, sorted_members)) + "}"
    print(result_str)

# Execute the function to find and print the result.
find_cycle_values()