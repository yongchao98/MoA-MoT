import sys

# It's good practice to increase the recursion limit for deep recursion,
# though it's not strictly necessary for n=10.
sys.setrecursionlimit(2000)

def solve():
    """
    Calculates the number of n-step self-avoiding walks on a Manhattan lattice.
    """
    n = 10

    def backtrack(steps_left, pos, visited):
        """
        Recursively counts the number of self-avoiding walks.

        Args:
            steps_left: The number of remaining steps in the walk.
            pos: The current (x, y) position on the lattice.
            visited: A set of (x, y) tuples representing visited points.

        Returns:
            The number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have found one valid walk.
        if steps_left == 0:
            return 1

        count = 0
        x, y = pos
        
        # Define the four possible moves from the current position.
        moves = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

        for next_pos in moves:
            # Check if the next position has already been visited.
            if next_pos not in visited:
                # If not visited, proceed with the walk from the new position.
                visited.add(next_pos)
                count += backtrack(steps_left - 1, next_pos, visited)
                # Backtrack: remove the point to explore other paths.
                visited.remove(next_pos)
        
        return count

    # The walk starts at the origin (0, 0).
    start_pos = (0, 0)
    # The initial set of visited points contains only the origin.
    visited = {start_pos}
    
    # Start the backtracking process.
    result = backtrack(n, start_pos, visited)
    
    # Print the final result in the format a(n) = result.
    print(f"a({n}) = {result}")

solve()