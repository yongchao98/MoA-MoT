def solve():
    """
    This script finds the total set of values included in fixed points or cycles
    for a specific digit reordering process, starting with all positive 
    three-digit numbers.
    """
    
    # A cache to memoize results of the calculation for efficiency.
    memo = {}

    def get_next_val(n):
        """
        Calculates the next value in the sequence for a given number n.
        The number is treated as having 3 digits, padding with leading zeros
        if its value is less than 100.
        """
        if n in memo:
            return memo[n]

        # Format n as a 3-digit string (e.g., 99 -> "099")
        s = f'{n:03d}'
        
        # Sort digits to find the smallest number A
        sorted_digits = sorted(list(s))
        a = int("".join(sorted_digits))
        
        # Reverse sorted digits to find the largest number B
        b = int("".join(sorted_digits[::-1]))
        
        result = b - a + 1
        memo[n] = result
        return result

    final_cycle_values = set()
    
    # Iterate through all positive three-digit numbers
    for start_num in range(100, 1000):
        # Track the path for the current sequence
        path = []
        visited_in_path = set()
        
        curr = start_num
        
        # Trace the sequence until a number is repeated, indicating a cycle.
        while curr not in visited_in_path:
            visited_in_path.add(curr)
            path.append(curr)
            curr = get_next_val(curr)
            
        # A cycle has been found. `curr` is the first repeated value.
        # Identify the numbers that are part of the cycle.
        try:
            cycle_start_index = path.index(curr)
            cycle = path[cycle_start_index:]
            # Add all numbers in the discovered cycle to our final set.
            final_cycle_values.update(cycle)
        except ValueError:
            # This case should not be reached due to the loop's logic.
            pass

    # Sort the final set of values in ascending order.
    sorted_values = sorted(list(final_cycle_values))
    
    # Print the result in the required format.
    print(f"{{{', '.join(map(str, sorted_values))}}}")

solve()