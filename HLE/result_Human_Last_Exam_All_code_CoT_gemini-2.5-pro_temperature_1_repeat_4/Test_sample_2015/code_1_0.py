def solve():
    """
    Solves the problem by simulating the process for all three-digit numbers,
    detecting cycles and fixed points, and collecting all values involved.
    """
    cycle_values = set()
    globally_visited = set()

    def process_step(n):
        """
        Calculates the next value in the sequence for a given number n.
        It finds the smallest number (A) and largest number (B) by reordering
        the digits of n, and returns B - A + 1.
        """
        s_n = str(n)
        
        # Form A, the smallest number
        a_digits = sorted(s_n)
        A = int("".join(a_digits))
        
        # Form B, the largest number
        b_digits = sorted(s_n, reverse=True)
        B = int("".join(b_digits))
        
        return B - A + 1

    # Iterate through all positive three-digit numbers
    for start_num in range(100, 1000):
        # If this number's path has been explored, skip it
        if start_num in globally_visited:
            continue

        path = []
        current_num = start_num

        # Trace the path until we find a cycle or hit a previously explored path
        while current_num not in globally_visited:
            # Check if we have entered a cycle within the current path
            if current_num in path:
                # A cycle is detected. Find where it starts.
                cycle_start_index = path.index(current_num)
                # The cycle consists of all elements from that index to the end
                the_cycle = path[cycle_start_index:]
                # Add all numbers in this newly found cycle to our result set
                cycle_values.update(the_cycle)
                break  # Stop tracing this path

            path.append(current_num)
            current_num = process_step(current_num)
        
        # Mark all numbers in the path we just traced as globally visited
        # to avoid re-computation in the future.
        globally_visited.update(path)

    # Sort the final set of values in ascending order
    sorted_cycle_values = sorted(list(cycle_values))
    
    # Format the output as a set string
    result_string = "{" + ", ".join(map(str, sorted_cycle_values)) + "}"
    print(result_string)

solve()
<<<1, 100, 397, 496, 595>>>