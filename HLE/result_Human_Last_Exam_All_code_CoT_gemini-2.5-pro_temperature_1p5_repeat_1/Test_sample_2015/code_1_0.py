import sys

def solve_cycles():
    """
    Finds the total set of values included in fixed values or cycles for a specific
    process run on all positive three-digit numbers.

    The process is:
    1. For a number n, reorder its digits to get the smallest (A) and largest (B) numbers.
    2. The next number is B - A + 1.
    3. Repeat this process.
    """

    # Memoization cache to store the cycle that a number eventually reaches.
    # Key: a number in a sequence, Value: the set of numbers in the terminal cycle.
    memo = {}

    # The set to store all unique values found in cycles or as fixed points.
    final_cycle_values = set()

    def get_next_value(n):
        """Calculates B - A + 1 for a given number n."""
        s = str(n)
        digits = sorted(list(s))
        
        # Form A, the smallest number
        str_a = "".join(digits)
        num_a = int(str_a)
        
        # Form B, the largest number
        str_b = "".join(reversed(digits))
        num_b = int(str_b)
        
        return num_b - num_a + 1

    # Iterate over all positive three-digit numbers (100 to 999)
    for start_num in range(100, 1000):
        # If we have already computed the cycle for this number, skip it.
        if start_num in memo:
            continue

        path = []
        current_num = start_num

        # Trace the sequence until we hit a number we've seen before,
        # either in the current path (new cycle) or in the memo (known cycle).
        while current_num not in memo and current_num not in path:
            path.append(current_num)
            current_num = get_next_value(current_num)

        # Determine the cycle.
        cycle_set = set()
        if current_num in memo:
            # We've reached a previously discovered cycle.
            cycle_set = memo[current_num]
        else:
            # A new cycle was found within the current path.
            cycle_start_index = path.index(current_num)
            cycle_members = path[cycle_start_index:]
            cycle_set = set(cycle_members)

        # Add the members of the found cycle to our final collection.
        final_cycle_values.update(cycle_set)

        # Memoize the result for every number in the path we just traced
        # to avoid re-computation.
        for num_in_path in path:
            memo[num_in_path] = cycle_set
            
    # Sort the final set of values in ascending order for the output.
    sorted_values = sorted(list(final_cycle_values))
    
    # Print the result in the specified format "{v1, v2, ...}".
    print(f"{{{', '.join(map(str, sorted_values))}}}")


solve_cycles()