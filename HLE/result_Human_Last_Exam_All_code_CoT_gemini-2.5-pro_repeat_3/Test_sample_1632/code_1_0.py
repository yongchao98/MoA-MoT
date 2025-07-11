def calculate_self_avoiding_walks(n_steps):
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice.

    This function uses a recursive backtracking algorithm to explore all possible
    self-avoiding paths of a given length starting from the origin (0,0).
    The set of visited coordinates is managed in the outer scope to be modified
    by the recursive helper function.
    """

    # The visited set stores the coordinates (x, y) of the path taken so far.
    # We start at the origin (0,0), so it's the first point in our path.
    visited = {(0, 0)}

    def _recursive_walk(current_n, x, y):
        """
        A helper function to perform the recursive search.
        
        Args:
            current_n: The number of steps remaining to take.
            x, y: The current coordinates on the lattice.
            
        Returns:
            The number of valid self-avoiding walks from the current position.
        """
        # Base case: If no steps are left, we have found one complete valid walk.
        if current_n == 0:
            return 1

        total_walks = 0
        
        # Explore the four possible directions: up, down, left, right.
        # The moves are represented as changes in x and y coordinates.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        
        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy
            
            # Check if the next position has already been visited.
            if (next_x, next_y) not in visited:
                # If not visited, explore this path.
                # 1. Add the new position to our current path.
                visited.add((next_x, next_y))
                
                # 2. Recursively call the function for the next step.
                total_walks += _recursive_walk(current_n - 1, next_x, next_y)
                
                # 3. Backtrack: Remove the new position so other paths can be explored.
                #    This is crucial for the algorithm to work correctly.
                visited.remove((next_x, next_y))
                
        return total_walks

    # Start the recursive search from the origin (0,0) with n_steps.
    return _recursive_walk(n_steps, 0, 0)

# The number of steps for the specific problem a(10).
n = 10

# Calculate the result by calling the main function.
result = calculate_self_avoiding_walks(n)

# Print the final result in the format "a(number) = result"
# as requested by "output each number in the final equation!".
print(f"a({n}) = {result}")