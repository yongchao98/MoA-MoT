def solve_self_avoiding_walks():
    """
    This script calculates a(n), the number of n-step self-avoiding walks
    on a 2D Manhattan (square) lattice, for n=10.
    """
    
    # The number of steps for the walk
    n = 10

    # A set to keep track of visited coordinates in the current path.
    # The walk starts at the origin, so it's the first point in our path.
    path = {(0, 0)}

    def count_walks(x, y, steps_left):
        """
        Recursively counts the number of valid self-avoiding walks starting
        from position (x, y) with a given number of steps remaining.
        
        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_left (int): The number of steps remaining in the walk.
            
        Returns:
            int: The total count of valid walks from the current state.
        """
        # Base case: If no steps are left, we have found one complete valid path.
        if steps_left == 0:
            return 1

        total_walks_from_here = 0
        
        # Explore the four neighbors on the Manhattan lattice (up, down, left, right)
        for dx, dy in [(0, 1), (0, -1), (-1, 0), (1, 0)]:
            next_x, next_y = x + dx, y + dy
            
            # Check if the next point is self-avoiding (not already in the path)
            if (next_x, next_y) not in path:
                # 1. Take a step: Add the new point to the path
                path.add((next_x, next_y))
                
                # 2. Recurse: Count walks from the new point with one less step
                total_walks_from_here += count_walks(next_x, next_y, steps_left - 1)
                
                # 3. Backtrack: Remove the point from the path so other recursive
                #    branches are not affected by this choice.
                path.remove((next_x, next_y))
                
        return total_walks_from_here

    # Initial call to start the recursive counting from the origin (0, 0)
    result = count_walks(0, 0, n)
    
    # Print the final result in the format a(n) = result
    print(f"a({n}) = {result}")

solve_self_avoiding_walks()