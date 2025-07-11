import sys

def solve():
    """
    Solves the problem of finding all cycle/fixed-point values for a given process
    run on all 3-digit numbers.
    """

    def get_next_val(n):
        """
        Calculates the next value in the sequence for a given number n.
        The next value is B - A + 1, where A and B are the smallest and largest
        numbers formed by the digits of n.
        """
        s_n = str(n)
        
        # Sort digits to find A
        sorted_digits = sorted(list(s_n))
        a = int("".join(sorted_digits))
        
        # Sort digits in reverse to find B
        # The list must be re-sorted from the original, not from the already sorted list
        reversed_digits = sorted(list(s_n), reverse=True)
        b = int("".join(reversed_digits))
        
        return b - a + 1

    cycle_nodes = set()
    # Memoization set to store nodes whose paths have already been fully explored.
    visited_nodes = set()

    # Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        # If this number's path has been explored, skip it.
        if start_num in visited_nodes:
            continue

        path = []
        path_seen = set()
        current_num = start_num

        # Trace the path until a cycle or an already visited node is found.
        while True:
            # If we hit a path that has been fully explored, we can stop.
            # We add the current path's nodes to the visited set for future optimization.
            if current_num in visited_nodes:
                visited_nodes.update(path)
                break

            # If we encounter a number already in the current path, we've found a cycle.
            if current_num in path_seen:
                try:
                    # The cycle starts from the first time we saw the current number.
                    cycle_start_index = path.index(current_num)
                    cycle = path[cycle_start_index:]
                    
                    # Add all numbers in the cycle to our final set.
                    cycle_nodes.update(cycle)
                    
                    # Add all nodes in the current path (tail + cycle) to the visited set.
                    visited_nodes.update(path)
                except ValueError:
                    # This should not be reached if current_num is in path_seen
                    pass
                break
            
            # Continue building the path.
            path.append(current_num)
            path_seen.add(current_num)
            current_num = get_next_val(current_num)

    # Sort the final set of cycle nodes for ascending order output.
    sorted_cycle_nodes = sorted(list(cycle_nodes))
    
    # Print the result in the specified format.
    result_str = "{" + ", ".join(map(str, sorted_cycle_nodes)) + "}"
    print(result_str)

solve()
<<<1, 100, 397, 496, 595>>>