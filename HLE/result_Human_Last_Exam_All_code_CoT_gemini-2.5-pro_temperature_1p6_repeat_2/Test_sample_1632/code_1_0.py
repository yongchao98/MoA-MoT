import sys

def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a 2D square lattice.
    """
    # Set a higher recursion limit for safety, as the path is 10 steps deep.
    sys.setrecursionlimit(2000)

    # The desired number of steps in the walk.
    N = 10

    def count_walks_recursive(x, y, steps_left, path):
        """
        Recursively counts self-avoiding walks using backtracking.

        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_left (int): The number of steps remaining in the walk.
            path (set): A set of (x, y) tuples representing visited points.

        Returns:
            int: The number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have successfully found one complete walk.
        if steps_left == 0:
            return 1

        count = 0
        # Explore the four possible directions (neighbors on a Manhattan lattice).
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_x, next_y = x + dx, y + dy

            # If the next point has not been visited, explore this path.
            if (next_x, next_y) not in path:
                # Add the new point to the path for the recursive call.
                path.add((next_x, next_y))
                # Recursively call the function for the remainder of the walk.
                count += count_walks_recursive(next_x, next_y, steps_left - 1, path)
                # Backtrack: remove the point to allow other paths to be explored.
                path.remove((next_x, next_y))
        return count

    # By symmetry, the total number of walks is 4 times the number of walks
    # that begin with a specific first step (e.g., moving right to (1, 0)).
    
    # We start the sub-problem from (1, 0), having already taken one step.
    # The path already contains the origin (0,0) and the first step (1,0).
    # The number of steps left to take is N - 1.
    initial_path = {(0, 0), (1, 0)}
    walks_starting_in_one_direction = count_walks_recursive(x=1, y=0, steps_left=N - 1, path=initial_path)
    
    # The total number of walks is 4 times this sub-problem result.
    total_walks = 4 * walks_starting_in_one_direction

    print(f"The number of {N}-step self-avoiding walks, a({N}), can be found using an optimized search.")
    print(f"By symmetry, we first calculate the number of walks that start by moving to (1, 0).")
    print(f"Number of walks starting with a single direction: {walks_starting_in_one_direction}")
    print("\nThe final equation is derived from multiplying this by the 4 possible starting directions:")
    print(f"a(10) = 4 * {walks_starting_in_one_direction} = {total_walks}")

solve()

<<<470092>>>