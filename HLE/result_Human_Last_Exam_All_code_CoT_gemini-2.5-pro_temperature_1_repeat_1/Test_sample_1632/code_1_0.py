def solve_a10():
    """
    Calculates a(10), the number of 10-step self-avoiding walks
    on a Manhattan lattice using a recursive backtracking algorithm.
    """
    n = 10
    
    # Memoization dictionary to store results of subproblems.
    # The key is a tuple: (steps_left, current_position, frozenset_of_relative_visited_points)
    # This optimization is complex to implement correctly for this problem structure.
    # A direct recursive search is feasible for n=10.
    
    def count_walks(steps_left, current_pos, visited):
        """
        Recursively counts the number of self-avoiding walks.

        Args:
            steps_left: The number of steps remaining in the walk.
            current_pos: The current (x, y) coordinate tuple.
            visited: A set of (x, y) tuples representing visited points.

        Returns:
            The total number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have found one complete valid path.
        if steps_left == 0:
            return 1

        x, y = current_pos
        total_walks = 0

        # Explore the four neighboring positions (North, South, East, West).
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_pos = (x + dx, y + dy)

            # If the next position has not been visited yet...
            if next_pos not in visited:
                # ...explore from this new position.
                visited.add(next_pos)
                total_walks += count_walks(steps_left - 1, next_pos, visited)
                # Backtrack: remove the point so other paths can be explored.
                visited.remove(next_pos)

        return total_walks

    # Initial call to the recursive function.
    # Start at the origin (0, 0) with n steps to go.
    # The initial path contains only the starting point.
    start_pos = (0, 0)
    visited_set = {start_pos}
    
    result = count_walks(n, start_pos, visited_set)

    # Print the final equation with the numbers involved.
    print(f"a({n}) = {result}")

# Execute the solution
solve_a10()
<<<44100>>>