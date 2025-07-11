def find_cycle_values():
    """
    This function calculates the set of all values that are part of a fixed point or cycle
    for a specific process applied to all 3-digit numbers.

    The process is as follows:
    1. For a number n, find A (smallest number from its digits) and B (largest).
    2. The next number is B - A + 1.
    3. Repeat with the new number.

    This function iterates through all starting numbers from 100 to 999, detects all
    cycles and fixed points, and collects their member values.
    """
    
    # Memoization cache to store the known cycle for a given number.
    memo = {}
    
    # A master set to store all numbers that are part of any cycle.
    master_cycle_set = set()

    def process_number(n):
        """
        Calculates the next number in the sequence from a given number n.
        Treats n as a 3-digit number, padding with leading zeros if needed.
        """
        # Pad with leading zeros to ensure 3 digits.
        s = str(n).zfill(3)
        
        # Sort digits to find the smallest (A) and largest (B) numbers.
        digits = sorted(list(s))
        a = int("".join(digits))
        b = int("".join(digits[::-1]))
        
        # Return the new value as per the problem description.
        return b - a + 1

    # Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        # If we have already processed this number's path, skip it.
        if start_num in memo:
            continue

        # Keep track of the current path to detect cycles.
        path = []
        # Use a set for faster lookups of numbers in the current path.
        path_set = set()
        current_num = start_num

        while True:
            # Case 1: The current number leads to a known cycle.
            if current_num in memo:
                cycle = memo[current_num]
                master_cycle_set.update(cycle)
                # Update memo for all numbers in the new path found.
                for p_num in path:
                    memo[p_num] = cycle
                break
            
            # Case 2: A new cycle is detected within the current path.
            if current_num in path_set:
                cycle_start_index = path.index(current_num)
                # The cycle consists of numbers from the first occurrence onwards.
                cycle_list = path[cycle_start_index:]
                cycle = set(cycle_list)
                
                master_cycle_set.update(cycle)
                
                # Update memo for all numbers in the path that lead to this cycle.
                for p_num in path:
                    memo[p_num] = cycle
                break
            
            # Continue building the path.
            path.append(current_num)
            path_set.add(current_num)
            current_num = process_number(current_num)

    # Sort the final set of values in ascending order.
    sorted_values = sorted(list(master_cycle_set))
    
    # Format the output as a set string.
    result_str = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_str)

find_cycle_values()