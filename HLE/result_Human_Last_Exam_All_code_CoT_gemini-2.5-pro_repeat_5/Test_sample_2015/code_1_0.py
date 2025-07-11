def solve_cycles():
    """
    This script finds the total set of values included in fixed points or cycles
    for a specific numerical process, starting with all positive three-digit numbers.
    The process is:
    1. For a number n, reorder its digits to get the smallest number A and largest number B.
    2. Calculate the new value as B - A + 1.
    3. Repeat the process with the new value.
    """

    def process_number(n):
        """
        Calculates the next number in the sequence (B - A + 1).
        The number of digits N is determined by the input n.
        """
        s_n = str(n)
        
        # Sort digits to get the smallest number, A.
        s_a = "".join(sorted(s_n))
        a = int(s_a)
        
        # Sort digits in reverse to get the largest number, B.
        s_b = "".join(sorted(s_n, reverse=True))
        b = int(s_b)
        
        return b - a + 1

    # This set will store all numbers that are part of a cycle or fixed point.
    final_set = set()
    
    # This set will keep track of numbers whose paths we've already explored.
    global_seen = set()

    # Iterate through all positive three-digit numbers.
    for start_n in range(100, 1000):
        # If we have already processed this number, skip it.
        if start_n in global_seen:
            continue
            
        path = []
        path_set = set()
        current_n = start_n
        
        # Trace the sequence until a number repeats.
        while current_n not in path_set:
            path.append(current_n)
            path_set.add(current_n)
            current_n = process_number(current_n)

        # Mark all numbers in the traced path as seen.
        global_seen.update(path)

        # A cycle has been found. `current_n` is the first repeated number.
        # We find where the cycle begins in our path.
        try:
            cycle_start_index = path.index(current_n)
            # The cycle consists of all numbers from the start of the repetition.
            cycle = path[cycle_start_index:]
            # Add all numbers from this cycle to our final set.
            final_set.update(cycle)
        except ValueError:
            # This block should not be reached due to the while loop's logic.
            pass

    # Sort the final list of numbers in ascending order for the output.
    sorted_values = sorted(list(final_set))
    
    # Print the result in the specified format, e.g., {1, 2, 3}
    print(f"{{{', '.join(map(str, sorted_values))}}}")

solve_cycles()