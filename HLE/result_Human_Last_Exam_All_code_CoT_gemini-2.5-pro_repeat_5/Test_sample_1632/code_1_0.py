def solve_self_avoiding_walks():
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice.
    This script is configured to calculate a(10).
    """
    n_steps = 10
    
    # A set to keep track of visited coordinates in the current path.
    # The walk starts at (0, 0), so we add it to the visited set.
    visited = {(0, 0)}

    def count_walks(x, y, steps_taken):
        """
        Recursively counts the number of self-avoiding walks using backtracking.
        
        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_taken (int): The number of steps taken so far.
        
        Returns:
            int: The number of valid self-avoiding walks from the current state.
        """
        # Base case: If we have taken the desired number of steps, we've found one valid walk.
        if steps_taken == n_steps:
            return 1

        count = 0
        
        # Explore the four possible directions (neighbors on a square lattice).
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_x, next_y = x + dx, y + dy
            
            # Check if the next point has already been visited.
            if (next_x, next_y) not in visited:
                # If not visited, explore this path.
                visited.add((next_x, next_y))
                count += count_walks(next_x, next_y, steps_taken + 1)
                
                # Backtrack: remove the point to explore other paths.
                visited.remove((next_x, next_y))
        
        return count

    # Start the search from the origin (0,0) with 0 steps taken.
    total_walks = count_walks(0, 0, 0)
    
    # The problem defines a(n) as the number of n-step walks.
    # The final output prints the result in the format a(n) = value.
    print(f"a({n_steps}) = {total_walks}")

solve_self_avoiding_walks()