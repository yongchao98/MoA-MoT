import sys

def solve_number_cycle():
    """
    Finds the total set of values included in fixed points or cycles for a specific
    number process applied to all 3-digit numbers.
    """

    def get_next_value(n):
        """Calculates B - A + 1 for a given number n."""
        # Pad with leading zeros to ensure a 3-digit string representation.
        s_n = str(n).zfill(3)
        
        # Sort digits to get the string for A.
        sorted_digits = sorted(s_n)
        s_A = "".join(sorted_digits)
        
        # Reverse sorted digits to get the string for B.
        s_B = "".join(sorted_digits[::-1])
        
        A = int(s_A)
        B = int(s_B)
        
        return B - A + 1

    # `visited` stores all numbers that have been part of any processed path.
    visited = set()
    # `cycle_values` will store the final set of numbers that are part of a cycle.
    cycle_values = set()

    # Iterate through all positive three-digit numbers as starting points.
    for i in range(100, 1000):
        # If we have already processed this number, skip it.
        if i in visited:
            continue
        
        # `path` stores the sequence generated from the current starting number.
        path = []
        current_num = i
        
        # Follow the path until we find a number that has been visited before.
        while current_num not in visited:
            visited.add(current_num)
            path.append(current_num)
            current_num = get_next_value(current_num)
            
        # At this point, current_num is the first number in the sequence that
        # was already visited. We need to check if it's part of the path
        # we just built, which would indicate a new cycle.
        if current_num in path:
            # A new cycle is found. Find where it starts in the current path.
            cycle_start_index = path.index(current_num)
            # Extract the cycle.
            new_cycle = path[cycle_start_index:]
            # Add all numbers from this new cycle to our final set.
            cycle_values.update(new_cycle)
            
    # Convert the set to a list and sort it in ascending order.
    sorted_cycle_values = sorted(list(cycle_values))
    
    # Print the result in the specified format "{v1, v2, ...}".
    result_str = f"{{{', '.join(str(v) for v in sorted_cycle_values)}}}"
    print(result_str)

solve_number_cycle()