def solve_number_cycle():
    """
    Finds the total set of values included in fixed points or cycles for a specific
    number transformation process applied to all 3-digit numbers.
    """

    def get_next_val(n):
        """
        Calculates the next value in the sequence for a given number n, assuming a 3-digit context.
        B is the largest number from n's digits, A is the smallest. Returns B - A + 1.
        """
        # Pad with leading zeros to ensure the context is always 3 digits.
        s_n = str(n).zfill(3)
        
        # Sort digits to get string representations for the smallest (s_A) and largest (s_B) numbers.
        digits = sorted(list(s_n))
        s_A = "".join(digits)
        s_B = "".join(reversed(digits))
        
        # Convert to integers to perform the calculation.
        A = int(s_A)
        B = int(s_B)
        
        return B - A + 1

    final_cycle_values = set()
    # A memoization cache to store the determined cycle for any number encountered.
    # This significantly speeds up the process by avoiding re-computation of paths.
    memo = {} 

    # Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        # If this number's path has already been computed, skip it.
        if start_num in memo:
            continue

        path = []
        current_num = start_num

        # Trace the path from the starting number.
        while current_num not in memo:
            # Check if the current number is already in the path we are building.
            # If so, we have found a new cycle.
            if current_num in path:
                cycle_start_index = path.index(current_num)
                cycle = path[cycle_start_index:]
                
                # Add the cycle's members to our final collection of values.
                final_cycle_values.update(cycle)
                
                # Memoize the result for all numbers on the path, as they all lead to this same cycle.
                canonical_cycle = tuple(sorted(cycle))
                for node in path:
                    memo[node] = canonical_cycle
                
                # Break the inner while loop since we've found and processed the cycle.
                break
            
            path.append(current_num)
            current_num = get_next_val(current_num)

        # This block handles the case where the path converges to a cycle that was already known.
        if start_num not in memo:
            # The current_num must be in the memo, so we've reached a known path.
            known_cycle = memo[current_num]
            # Memoize all nodes in the newly traced path to point to this known cycle.
            for node in path:
                memo[node] = known_cycle
                
    # Sort the final set of values in ascending order.
    sorted_values = sorted(list(final_cycle_values))
    
    # Format the output as a set string, e.g., {1, 2, 3}.
    output_str = f"{{{', '.join(map(str, sorted_values))}}}"
    print(output_str)

solve_number_cycle()