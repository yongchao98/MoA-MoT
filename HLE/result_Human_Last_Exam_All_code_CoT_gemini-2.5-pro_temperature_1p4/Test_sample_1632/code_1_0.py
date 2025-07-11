def calculate_self_avoiding_walks():
    """
    Calculates the number of n-step self-avoiding walks on a Manhattan lattice
    using a recursive backtracking algorithm.
    """
    # The desired length of the walk.
    n = 10
    
    # A counter for the number of valid walks found.
    # We use a list to make it mutable inside the nested function.
    walk_count = [0]

    def explore(current_pos, path, length):
        """
        Recursively explores paths from the current position.

        Args:
            current_pos (tuple): The (x, y) coordinates of the current step.
            path (set): A set of (x, y) tuples representing visited points.
            length (int): The current length of the walk.
        """
        # Base case: If the walk has reached the desired length, increment the counter.
        if length == n:
            walk_count[0] += 1
            return

        x, y = current_pos
        
        # Define the four possible moves: North, East, South, West.
        moves = [(x, y + 1), (x + 1, y), (x, y - 1), (x - 1, y)]

        # Explore each possible move.
        for next_pos in moves:
            # A move is valid if the next position has not been visited.
            if next_pos not in path:
                # Add the new position to the path and recurse.
                path.add(next_pos)
                explore(next_pos, path, length + 1)
                # Backtrack: remove the position to explore other branches.
                path.remove(next_pos)

    # The walk starts at the origin (0, 0).
    start_position = (0, 0)
    
    # The initial path contains only the starting point.
    initial_path = {start_position}

    # Start the recursive exploration.
    explore(start_position, initial_path, 0)
    
    # The final result is the first element of our counter list.
    result = walk_count[0]
    
    # Print the result, including the numbers from the equation a(n) = result.
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    calculate_self_avoiding_walks()