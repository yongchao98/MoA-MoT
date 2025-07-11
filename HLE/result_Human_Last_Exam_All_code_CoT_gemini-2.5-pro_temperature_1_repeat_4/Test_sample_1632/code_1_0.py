def solve_saw_count():
    """
    Calculates the number of n-step self-avoiding walks on a 2D square lattice.
    """
    n = 10

    def count_walks_recursive(x, y, steps_left, visited):
        """
        Recursively counts the number of self-avoiding walks.

        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_left (int): The number of steps remaining in the walk.
            visited (set): A set of (x, y) tuples representing visited points.

        Returns:
            int: The number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have successfully found one valid walk.
        if steps_left == 0:
            return 1

        count = 0
        # Explore the four neighboring points.
        # Moves are represented as (dx, dy).
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nx, ny = x + dx, y + dy

            # Check if the neighbor has been visited.
            if (nx, ny) not in visited:
                # If not visited, add it to the path and recurse.
                visited.add((nx, ny))
                count += count_walks_recursive(nx, ny, steps_left - 1, visited)
                # Backtrack: remove the point to explore other paths.
                visited.remove((nx, ny))
        
        return count

    # The walk starts at the origin (0, 0).
    start_pos = (0, 0)
    # The initial path contains only the starting point.
    initial_visited = {start_pos}

    # Start the recursive counting process.
    total_walks = count_walks_recursive(start_pos[0], start_pos[1], n, initial_visited)

    # Print the final result in the format a(n) = result.
    print(f"a({n}) = {total_walks}")

solve_saw_count()