def solve_kaprekar_cycle():
    """
    This script finds the set of all numbers that are part of a fixed point or cycle
    for a specific mathematical process.

    The process is as follows:
    1. Take a number n.
    2. Form the largest number B and smallest number A from its digits.
    3. The next number is B - A + 1.
    4. Repeat this process.

    This script runs the process for all positive three-digit numbers (100-999)
    and identifies the complete set of values that are part of any resulting
    cycle or fixed point.
    """
    
    # memo will store the terminal cycle for any number encountered.
    # e.g., memo[n] = {cycle_val_1, cycle_val_2, ...}
    memo = {}
    
    def get_next_value(n):
        """Calculates B - A + 1 for a given number n."""
        s_n = str(n)
        s_a = "".join(sorted(s_n))
        s_b = "".join(sorted(s_n, reverse=True))
        a = int(s_a)
        b = int(s_b)
        return b - a + 1

    def find_cycle_for_start_node(start_node):
        """
        Traces the path from a starting number until a cycle is found.
        Uses and populates a global memo to avoid re-computation.
        """
        # If the cycle for this node is already known, return it.
        if start_node in memo:
            return memo[start_node]

        path_sequence = []
        path_set = set()
        current_node = start_node

        # Traverse the path until we find a node that is already in our
        # current path or has a known cycle in the memo.
        while current_node not in path_set and current_node not in memo:
            path_set.add(current_node)
            path_sequence.append(current_node)
            current_node = get_next_value(current_node)

        # Determine the cycle
        if current_node in memo:
            # The path leads to a previously computed cycle
            cycle = memo[current_node]
        else:
            # A new cycle is found within the current path
            cycle_start_index = path_sequence.index(current_node)
            cycle = set(path_sequence[cycle_start_index:])
        
        # Memoize the found cycle for all nodes in the current path
        for node in path_sequence:
            memo[node] = cycle
            
        return cycle

    # A set to store all unique numbers found in any cycle
    all_cycle_values = set()

    # Iterate through all three-digit numbers
    for i in range(100, 1000):
        cycle = find_cycle_for_start_node(i)
        all_cycle_values.update(cycle)
        
    # Sort the final set of values and print in the required format
    sorted_values = sorted(list(all_cycle_values))
    
    # The string formatting below creates the output like {1, 2, 3}
    result_string = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_string)

solve_kaprekar_cycle()
<<<1, 100, 397, 496, 595>>>