def count_self_avoiding_walks(n):
    """
    Calculates the number of n-step self-avoiding walks on a 2D lattice.
    
    This function uses a backtracking algorithm to explore all possible paths
    starting from the origin (0,0).
    """
    # A set to store the coordinates of visited points for efficient lookup.
    # The walk starts at (0,0), so it's the first visited point.
    visited = set([(0, 0)])

    def backtrack(x, y, steps_remaining):
        """
        Recursively counts the number of valid walks from the current position.
        
        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_remaining (int): The number of steps left to take.
            
        Returns:
            int: The number of valid self-avoiding walks from this point.
        """
        # Base case: If there are no steps left, we have found one complete walk.
        if steps_remaining == 0:
            return 1
        
        # Initialize a counter for the number of valid walks from this point.
        walk_count = 0
        
        # Define the four possible moves on a Manhattan lattice (up, down, left, right).
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        
        # Explore each possible move.
        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy
            
            # Check if the next point has already been visited.
            if (next_x, next_y) not in visited:
                # 'Choose' this move: add the new point to the visited set.
                visited.add((next_x, next_y))
                
                # 'Explore': Recursively call the function for the next step.
                walk_count += backtrack(next_x, next_y, steps_remaining - 1)
                
                # 'Unchoose' (backtrack): Remove the point to allow other paths to use it.
                visited.remove((next_x, next_y))
                
        return walk_count

    # Start the backtracking process from the origin (0,0) with n steps.
    return backtrack(0, 0, n)

# The number of steps for the walk.
n = 10

# Calculate the result.
result = count_self_avoiding_walks(n)

# Print the final result in the format "a(n) = result".
print(f"a({n}) = {result}")
