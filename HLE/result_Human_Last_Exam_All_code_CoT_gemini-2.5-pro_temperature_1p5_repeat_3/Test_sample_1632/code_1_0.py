def solve_saw_count():
    """
    Calculates and prints the number of 10-step self-avoiding walks 
    on a Manhattan lattice.
    """
    n = 10

    def count_walks(x, y, steps, visited):
        """
        Recursively counts valid walks from position (x,y) using backtracking.
        
        Args:
            x (int): current x-coordinate.
            y (int): current y-coordinate.
            steps (int): number of steps taken so far.
            visited (set): a set of (x, y) tuples representing the current path.

        Returns:
            int: The number of valid self-avoiding walks from this state.
        """
        # Base case: If we have completed n steps, we've found one valid walk.
        if steps == n:
            return 1

        count = 0
        # The four possible moves on a Manhattan lattice: East, West, North, South.
        moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        
        # Explore each possible move.
        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy
            
            # Check if the next position has already been visited.
            if (next_x, next_y) not in visited:
                # If not, add the new position to the path.
                visited.add((next_x, next_y))
                # Recursively call for the next step.
                count += count_walks(next_x, next_y, steps + 1, visited)
                # Backtrack: remove the position to explore other paths.
                visited.remove((next_x, next_y))
        
        return count

    # The walk starts at the origin (0, 0).
    start_x, start_y = 0, 0
    # The initial path contains only the starting point.
    initial_visited = {(start_x, start_y)}
    
    # Start the calculation.
    result = count_walks(start_x, start_y, 0, initial_visited)
    
    # Print the result in the requested equation format.
    print(f"a({n}) = {result}")

solve_saw_count()