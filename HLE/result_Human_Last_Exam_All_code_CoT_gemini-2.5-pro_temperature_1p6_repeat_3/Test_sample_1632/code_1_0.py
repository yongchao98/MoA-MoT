def solve_a10():
    """
    Calculates and prints a(10), the number of 10-step self-avoiding walks
    on a Manhattan (2D square) lattice.
    """
    
    n = 10  # The number of steps in the walk

    def count_walks(steps, path, x, y):
        """
        Recursively counts the number of self-avoiding walks.

        Args:
            steps (int): The number of steps remaining.
            path (set): A set of (x, y) tuples representing visited coordinates.
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.

        Returns:
            int: The number of valid walks from the current state.
        """
        # Base case: If no steps are left, we have found one valid walk.
        if steps == 0:
            return 1

        count = 0
        # Define the four possible moves: up, down, left, right.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        # Explore each possible move.
        for dx, dy in moves:
            nx, ny = x + dx, y + dy

            # Check if the next position has already been visited.
            if (nx, ny) not in path:
                # If not visited, take the step.
                path.add((nx, ny))
                
                # Recursively call the function for the path with one less step.
                count += count_walks(steps - 1, path, nx, ny)
                
                # Backtrack: remove the point from the path to explore other possibilities.
                path.remove((nx, ny))
        
        return count

    # The walk starts at the origin (0, 0).
    start_x, start_y = 0, 0
    
    # The initial path contains only the origin.
    # A set is used for efficient checking of visited points (O(1) average time complexity).
    initial_path = {(start_x, start_y)}

    # Calculate the total number of walks.
    result = count_walks(n, initial_path, start_x, start_y)
    
    # The final equation is a(10) = result.
    # We print each number involved in this equation.
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    solve_a10()