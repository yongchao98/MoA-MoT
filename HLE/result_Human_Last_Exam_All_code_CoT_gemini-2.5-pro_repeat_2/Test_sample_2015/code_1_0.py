def find_all_cycle_values():
    """
    Finds the total set of values in fixed points or cycles for the given process
    run on all positive three-digit numbers.
    """

    def get_next_value(n):
        """
        Calculates the next value in the sequence for a given number n.
        The number is treated as a 3-digit number.
        """
        # Pad with leading zeros to ensure 3 digits
        s_n = str(n).zfill(3)
        
        # Sort digits to get components for A and B
        digits = sorted(list(s_n))
        
        # Form smallest number A and largest number B
        s_A = "".join(digits)
        s_B = "".join(reversed(digits))
        
        A = int(s_A)
        B = int(s_B)
        
        return B - A + 1

    # This set will store all numbers that are part of a cycle or are fixed points.
    all_cycle_members = set()
    
    # Memoization cache to store the cycle associated with a number.
    # memo[number] = {cycle_members}
    memo = {}

    # Iterate through all positive three-digit numbers.
    for i in range(100, 1000):
        # Skip if this number's path has already been computed.
        if i in memo:
            continue

        path = []
        current_num = i

        # Trace the path until a cycle is found or a known path is hit.
        while True:
            if current_num in memo:
                # We've hit a number whose cycle is already known.
                cycle = memo[current_num]
                all_cycle_members.update(cycle)
                # Memoize the path that led here.
                for node in path:
                    if node not in memo:
                        memo[node] = cycle
                break

            if current_num in path:
                # We've discovered a new cycle.
                cycle_start_index = path.index(current_num)
                cycle_nodes = path[cycle_start_index:]
                cycle = set(cycle_nodes)
                
                all_cycle_members.update(cycle)
                
                # Memoize for all nodes in the path leading to and including the cycle.
                for node in path:
                    if node not in memo:
                        memo[node] = cycle
                break
            
            path.append(current_num)
            current_num = get_next_value(current_num)

    # Sort the final set for the required output format.
    sorted_values = sorted(list(all_cycle_members))
    
    # Format the output string as {v1, v2, v3, ...}
    result_string = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_string)

find_all_cycle_values()
<<< {100, 397, 496, 595} >>>