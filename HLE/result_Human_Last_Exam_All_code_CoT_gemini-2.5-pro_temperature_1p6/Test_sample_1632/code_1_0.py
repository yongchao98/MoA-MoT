def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice.
    """
    
    # The number of steps for the walk.
    n = 10

    def count_walks(steps, x, y, path):
        """
        Recursively counts the number of self-avoiding walks using backtracking.
        
        Args:
            steps: The number of steps remaining.
            x, y: The current coordinates.
            path: A set of (x, y) tuples of visited points.
            
        Returns:
            The number of valid walks from the current state.
        """
        # Base case: If no steps are left, we have found one complete walk.
        if steps == 0:
            return 1

        total_walks = 0
        
        # The four possible moves on a square lattice.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        for dx, dy in moves:
            new_x, new_y = x + dx, y + dy
            
            # Check if the next point is already in the path.
            if (new_x, new_y) not in path:
                # If not visited, add it to the path and recurse.
                path.add((new_x, new_y))
                total_walks += count_walks(steps - 1, new_x, new_y, path)
                
                # Backtrack: remove the point to explore other paths.
                path.remove((new_x, new_y))
        
        return total_walks

    # The walk starts at the origin (0, 0).
    start_x, start_y = 0, 0
    
    # The initial path contains only the starting point.
    initial_path = {(start_x, start_y)}

    # Start the calculation.
    result = count_walks(n, start_x, start_y, initial_path)

    # Print the final result as an equation.
    print(f"a({n}) = {result}")

solve()