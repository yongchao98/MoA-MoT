def solve_saw():
    """
    Calculates the number of n-step self-avoiding walks on a 2D Manhattan lattice.
    """
    # The problem is to find a(10).
    n = 10

    def count_walks(steps_left, x, y, visited):
        """
        Recursively counts self-avoiding walks using backtracking.

        Args:
            steps_left (int): The number of steps remaining.
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            visited (set): A set of (x, y) tuples of visited points.

        Returns:
            int: The number of valid walks from the current state.
        """
        # Base case: If no steps are left, we've found one valid walk.
        if steps_left == 0:
            return 1

        count = 0
        # The four possible moves: North, South, East, West.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy

            # A move is valid only if the next point has not been visited.
            if (next_x, next_y) not in visited:
                # Mark the new point as visited for this path.
                visited.add((next_x, next_y))
                # Recursively count walks from the new point.
                count += count_walks(steps_left - 1, next_x, next_y, visited)
                # Backtrack: remove the point to allow other paths to be explored.
                visited.remove((next_x, next_y))

        return count

    # The walk starts at the origin (0, 0).
    start_point = (0, 0)
    # The initial set of visited points contains only the origin.
    initial_visited = {start_point}

    # Start the calculation.
    total_walks = count_walks(n, start_point[0], start_point[1], initial_visited)

    # As requested, output the numbers in the final equation a(n) = result.
    print(f"a({n}) = {total_walks}")

solve_saw()