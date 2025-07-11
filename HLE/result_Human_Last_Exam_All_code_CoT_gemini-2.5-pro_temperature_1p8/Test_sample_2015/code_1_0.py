import collections

def find_cycle_values():
    """
    Runs the specified process on all positive three-digit numbers to find the total
    set of values included in the final fixed values or cycles.
    """
    # This set will store all unique numbers that are part of any cycle.
    final_cycle_values = set()

    def get_next_value(n):
        """Calculates the next value in the sequence for a given number n."""
        # The problem is defined for N=3 digits, so we pad with zeros if needed.
        s_n = str(n).zfill(3)
        
        # Sort digits to find A (smallest) and B (largest).
        sorted_digits = sorted(list(s_n))
        
        str_a = "".join(sorted_digits)
        str_b = "".join(reversed(sorted_digits))
        
        a = int(str_a)
        b = int(str_b)
        
        return b - a + 1

    # Iterate through every positive three-digit number.
    for start_n in range(100, 1000):
        
        current_n = start_n
        history_list = []
        history_set = set()

        # Follow the sequence until a number repeats.
        while current_n not in history_set:
            history_list.append(current_n)
            history_set.add(current_n)
            current_n = get_next_value(current_n)

        # A cycle is found. The repeated number is 'current_n'.
        # Find where the cycle begins in our history list.
        try:
            cycle_start_index = history_list.index(current_n)
            # Extract the members of the cycle.
            cycle = history_list[cycle_start_index:]
            # Add all numbers from this cycle to our final collection.
            final_cycle_values.update(cycle)
        except ValueError:
            # This case should not be reached with the current logic.
            pass

    # Sort the final results for the required output format.
    sorted_values = sorted(list(final_cycle_values))
    
    # Print the set in the specified format.
    # The f-string automatically converts list elements to strings.
    print(f"{{{', '.join(map(str, sorted_values))}}}")


find_cycle_values()