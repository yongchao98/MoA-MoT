import sys

def solve_saw_count():
    """
    Calculates and prints the number of n-step self-avoiding walks on a 2D lattice.
    """
    # The number of steps for the walk.
    N = 10

    # It's possible for the recursion to go deep. We can increase the Python's
    # recursion limit, although the default is usually sufficient for N=10.
    # sys.setrecursionlimit(2000)

    # Use a set for the visited path for efficient O(1) lookups.
    visited_path = set()

    def count_walks(x, y, steps_left):
        """
        Recursively counts the number of self-avoiding walks using backtracking.

        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_left (int): The number of steps remaining in the walk.

        Returns:
            int: The number of valid self-avoiding walks from the current state.
        """
        # Base case: If there are no steps left, we have found one valid path.
        if steps_left == 0:
            return 1

        # Add the current point to the visited path.
        visited_path.add((x, y))

        # Initialize the count of walks from this point.
        walk_count = 0

        # Define the four possible moves on a Manhattan lattice.
        moves = [(0, 1),  # Up
                 (0, -1), # Down
                 (1, 0),  # Right
                 (-1, 0)] # Left

        # Explore each possible move.
        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy

            # If the next point has not been visited, explore from there.
            if (next_x, next_y) not in visited_path:
                walk_count += count_walks(next_x, next_y, steps_left - 1)

        # Backtrack: Remove the current point from the path so it can be
        # visited by other branches of the search.
        visited_path.remove((x, y))

        return walk_count

    # The walk starts at the origin (0,0).
    start_x, start_y = 0, 0

    # The initial call to the function to count all walks of length N.
    total_walks = count_walks(start_x, start_y, N)

    # The provided prompt specifies an unusual output format for the "final equation"
    # The clearest way to present this is to show the value of a(10).
    print(f"a({N}) = {total_walks}")


solve_saw_count()