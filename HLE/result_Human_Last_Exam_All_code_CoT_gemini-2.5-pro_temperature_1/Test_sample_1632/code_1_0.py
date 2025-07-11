def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks
    on a Manhattan lattice.
    """
    N = 10

    def count_walks_recursive(steps_left, x, y, visited):
        """
        Recursively counts the number of self-avoiding walks from a given state.

        Args:
            steps_left: The number of steps remaining in the walk.
            x, y: The coordinates of the current position.
            visited: A set of (x, y) tuples representing visited points.

        Returns:
            The number of valid walks from the current state.
        """
        # Base case: If no steps are left, we have found one complete walk.
        if steps_left == 0:
            return 1

        count = 0
        # The four possible moves on a Manhattan lattice.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy
            # Check if the next point has already been visited.
            if (next_x, next_y) not in visited:
                # If not visited, add it to the path and explore further.
                visited.add((next_x, next_y))
                count += count_walks_recursive(steps_left - 1, next_x, next_y, visited)
                # Backtrack: remove the point so other branches of the search can use it.
                visited.remove((next_x, next_y))
        
        return count

    # Main calculation starts here.
    # By symmetry, we can calculate walks for one initial direction and multiply by 4.
    # We choose the first step to be from (0,0) to (1,0).
    # This means we have N-1 steps left to calculate.
    
    start_x, start_y = 0, 0
    first_step_x, first_step_y = 1, 0
    
    # The visited set starts with the origin and the point after the first step.
    initial_visited = {(start_x, start_y), (first_step_x, first_step_y)}
    
    # Calculate the number of walks for one fixed initial direction.
    walks_one_direction = count_walks_recursive(N - 1, first_step_x, first_step_y, initial_visited)
    
    # The total number of walks is 4 times this value.
    a_10 = 4 * walks_one_direction

    print(f"Let a(n) be the number of n-step self-avoiding walks on a Manhattan lattice.")
    print(f"To find a(10), we count all valid paths of length 10 starting from the origin.")
    print(f"Using symmetry, we calculate paths starting with one specific move (e.g., right) and multiply by 4.")
    print(f"Number of paths starting with a step to the right = {walks_one_direction}")
    print(f"a(10) = 4 * {walks_one_direction}")
    print(f"a(10) = {a_10}")

solve()
<<<22144>>>