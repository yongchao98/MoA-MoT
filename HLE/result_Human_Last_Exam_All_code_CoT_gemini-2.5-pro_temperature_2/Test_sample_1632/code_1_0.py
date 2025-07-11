import sys

def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice.
    """
    
    # On some systems, the default recursion limit might be too low for deep searches,
    # though it should be sufficient for n=10.
    # It's good practice to increase it for such problems.
    sys.setrecursionlimit(2000)

    n = 10

    # The recursive backtracking function to count walks.
    # A dictionary is used as a simple cache (memoization), though for n=10, 
    # the performance gain is not strictly necessary but good practice.
    # The state is defined by the remaining steps and a frozenset of the path
    # relative to the current position. This is complex to implement correctly,
    # so we will use a simpler direct recursion.
    def count_walks(steps_left, x, y, visited):
        """
        Recursively counts the number of self-avoiding walks.

        Args:
            steps_left (int): Number of steps remaining.
            x (int): Current x-coordinate.
            y (int): Current y-coordinate.
            visited (set): A set of (x, y) tuples of visited points.

        Returns:
            int: The number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have found one valid path.
        if steps_left == 0:
            return 1

        total_walks = 0
        
        # Explore the four possible directions (Manhattan lattice).
        moves = [(0, 1),  # Up
                 (0, -1), # Down
                 (1, 0),  # Right
                 (-1, 0)] # Left

        for dx, dy in moves:
            new_x, new_y = x + dx, y + dy

            # Check if the new point has been visited.
            if (new_x, new_y) not in visited:
                # If not, explore from this new point.
                visited.add((new_x, new_y))
                total_walks += count_walks(steps_left - 1, new_x, new_y, visited)
                
                # Backtrack: remove the point to explore other paths.
                visited.remove((new_x, new_y))

        return total_walks

    # Initial setup for the walk
    start_pos = (0, 0)
    # The initial path contains only the starting point.
    initial_visited = {start_pos}

    # Start the calculation
    result = count_walks(n, start_pos[0], start_pos[1], initial_visited)
    
    # Print the final equation with the numbers and result.
    print(f"Let a(n) be the number of n-step self-avoiding walks on a Manhattan lattice.")
    print(f"To find a(10), we count all valid paths of 10 steps.")
    print(f"a({n}) = {result}")

solve()
<<<44100>>>