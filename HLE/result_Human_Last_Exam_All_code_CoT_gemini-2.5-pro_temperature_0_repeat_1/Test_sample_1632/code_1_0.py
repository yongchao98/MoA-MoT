def solve_self_avoiding_walks():
    """
    Calculates a(n), the number of n-step self-avoiding walks on a 
    Manhattan (2D square) lattice.
    """
    n = 10

    def count_walks_recursive(steps_left, current_pos, visited):
        """
        Recursively counts the number of self-avoiding walks from a given state.

        Args:
            steps_left: The number of steps remaining in the walk.
            current_pos: The current (x, y) coordinates on the lattice.
            visited: A set of (x, y) tuples representing visited points.

        Returns:
            The number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have successfully found one complete walk.
        if steps_left == 0:
            return 1

        count = 0
        x, y = current_pos

        # Explore the four neighbors: up, down, right, left.
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_pos = (x + dx, y + dy)

            # If the neighbor has not been visited, proceed with the walk.
            if next_pos not in visited:
                # Mark the neighbor as visited for the subsequent recursive calls.
                visited.add(next_pos)
                # Recursively call the function for the next step.
                count += count_walks_recursive(steps_left - 1, next_pos, visited)
                # Backtrack: Un-mark the neighbor as visited to explore other paths.
                visited.remove(next_pos)
        
        return count

    # A walk of 0 steps is just the starting point, so there's 1 such walk.
    if n == 0:
        result = 1
        print(f"a(0) = {result}")
        return

    # Start the walk at the origin (0, 0).
    start_pos = (0, 0)
    
    # Optimization: Exploit the symmetry of the grid. We calculate the number of
    # walks starting with a step in one direction (e.g., right) and multiply by 4.
    # The first step is from (0, 0) to (1, 0).
    first_step_pos = (1, 0)
    
    # The set of visited points initially contains the start and the first step.
    initial_visited = {start_pos, first_step_pos}
    
    # We now need to find the number of walks with n-1 steps remaining.
    walks_in_one_direction = count_walks_recursive(n - 1, first_step_pos, initial_visited)
    
    # The total number of walks is 4 times this amount.
    total_walks = 4 * walks_in_one_direction

    print(f"Let a(n) be the number of n-step self-avoiding walks.")
    print(f"We want to find a({n}).")
    print(f"By symmetry, we can calculate the number of walks starting with a step to the right and multiply by 4.")
    print(f"Number of walks starting with one step to the right: {walks_in_one_direction}")
    print(f"a({n}) = 4 * {walks_in_one_direction} = {total_walks}")

# Execute the solution
solve_self_avoiding_walks()