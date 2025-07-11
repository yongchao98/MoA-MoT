def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice,
    using a recursive backtracking algorithm with a symmetry optimization.
    """

    # The total number of steps for the walk.
    total_steps = 10

    def count_walks_recursive(steps_remaining, current_pos, path):
        """
        Recursively counts valid self-avoiding walks from a given state.

        Args:
            steps_remaining: Number of steps left to take.
            current_pos: The current (x, y) coordinates of the walk.
            path: A set of (x, y) tuples representing the visited points.

        Returns:
            The total number of valid self-avoiding walks from the current state.
        """
        # Base case: If there are no steps left, we have found one valid walk.
        if steps_remaining == 0:
            return 1

        x, y = current_pos
        count = 0

        # Explore the 4 possible moves: right, left, up, down.
        for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            next_pos = (x + dx, y + dy)
            
            # Check if the next position has already been visited.
            if next_pos not in path:
                # If the position is new, continue the walk from there.
                # We pass a new set `path | {next_pos}` to the recursive call.
                count += count_walks_recursive(
                    steps_remaining - 1, 
                    next_pos, 
                    path | {next_pos}
                )
        
        return count

    # Optimization: Due to lattice symmetry, we can calculate the number of walks
    # for one initial direction (e.g., right) and multiply the result by 4.
    # The walk starts at (0,0), and the first step is to (1,0).
    # We then need to count the number of (total_steps - 1) walks from (1,0).
    
    # The starting point of the walk.
    origin = (0, 0)
    # The position after the first fixed step to the right.
    first_step_pos = (1, 0)
    
    # The initial path for the recursion already contains the origin and the first step.
    initial_path = {origin, first_step_pos}
    
    # Number of steps remaining after the first one.
    steps_after_first = total_steps - 1

    # Calculate the number of walks for one starting direction.
    count_one_direction = count_walks_recursive(
        steps_after_first, 
        first_step_pos, 
        initial_path
    )
    
    # The total number of walks is 4 times the count for one direction.
    total_walks = 4 * count_one_direction
    
    print(f"To find a({total_steps}), we calculate the number of walks for one initial direction and multiply by 4.")
    print(f"Number of {total_steps}-step walks starting with a move to the right: {count_one_direction}")
    print(f"a({total_steps}) = 4 * {count_one_direction} = {total_walks}")

# Execute the solver function
solve()
<<<44100>>>