import collections

def find_cycle_values():
    """
    Finds all numbers that are part of a fixed point or cycle for a given process
    on three-digit numbers.
    """
    final_cycle_values = set()

    def process_step(n):
        """
        Performs the operation: n -> B - A + 1, where A and B are the smallest
        and largest numbers formed by the digits of n.
        """
        # Ensure n is treated as a 3-digit number by padding with zeros.
        s = str(n).zfill(3)
        
        # Create the smallest number A.
        a = int("".join(sorted(s)))
        
        # Create the largest number B.
        b = int("".join(sorted(s, reverse=True)))
        
        return b - a + 1

    # Iterate over all positive three-digit numbers.
    for start_num in range(100, 1000):
        
        path = []
        visited_in_path = set()
        current_num = start_num
        
        # Follow the sequence until a number repeats in the current path.
        while current_num not in visited_in_path:
            visited_in_path.add(current_num)
            path.append(current_num)
            current_num = process_step(current_num)
            
        # A cycle has been found. The cycle begins at the first occurrence
        # of the repeated number `current_num` in the path list.
        try:
            cycle_start_index = path.index(current_num)
            cycle = path[cycle_start_index:]
            final_cycle_values.update(cycle)
        except ValueError:
            # This should not be reached with the current logic.
            pass

    # Sort the final set of values and print in the specified format.
    sorted_results = sorted(list(final_cycle_values))
    
    # Format the output string as requested, e.g., {1, 2, 3}
    output_str = "{" + ", ".join(map(str, sorted_results)) + "}"
    print(output_str)

find_cycle_values()