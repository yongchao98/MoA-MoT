def solve_saw_problem():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice.

    A self-avoiding walk is a path on a lattice that does not visit the same point
    more than once. This function uses a recursive backtracking algorithm to count
    all such walks of a specified length.

    To optimize, it uses the symmetry of the square lattice. The total number of
    walks is 4 times the number of walks starting in a specific direction (e.g., East).
    """
    n = 10

    def count_walks_recursive(path, x, y, length):
        """
        Recursively counts self-avoiding walks from a given state.

        Args:
            path: A set of (x, y) tuples representing visited points.
            x, y: The current coordinates of the head of the walk.
            length: The current number of steps taken.

        Returns:
            The number of valid n-step self-avoiding walks from the current state.
        """
        # Base case: If the walk has reached the desired length, we found one valid walk.
        if length == n:
            return 1

        count = 0
        # Explore the four possible moves (North, South, East, West).
        moves = [(x, y + 1), (x, y - 1), (x + 1, y), (x - 1, y)]

        for next_x, next_y in moves:
            # Check if the next point has already been visited.
            if (next_x, next_y) not in path:
                # If not, add the new point to the path and recurse.
                path.add((next_x, next_y))
                count += count_walks_recursive(path, next_x, next_y, length + 1)
                # Backtrack: remove the point to allow exploration of other branches.
                path.remove((next_x, next_y))
        return count

    # Start the walk at (0,0).
    # To use symmetry, we fix the first step, for example, to (1,0).
    # The path now contains two points, and the length of the walk is 1.
    start_path = {(0, 0), (1, 0)}
    
    # We need to find the number of walks of length n-1 starting from (1,0),
    # with (0,0) already visited.
    count_one_direction = count_walks_recursive(start_path, 1, 0, 1)
    
    # The total number of walks is 4 times this count.
    total_walks = 4 * count_one_direction

    print(f"For n = {n}:")
    print(f"Number of walks starting in one fixed direction = {count_one_direction}")
    print(f"a({n}) = 4 * {count_one_direction} = {total_walks}")

solve_saw_problem()
<<<44100>>>