import sys

def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice.
    """
    N = 10
    
    # It's good practice to have a higher recursion limit for deep recursion,
    # though for N=10, the default limit is typically sufficient.
    # We increase it to be safe.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), N + 5))

    def count_walks_recursive(steps_left, current_pos, visited):
        """
        Recursively counts the number of valid walk extensions from a given state.
        
        Args:
            steps_left: Number of steps remaining in the walk.
            current_pos: The current (x, y) coordinate tuple.
            visited: A set of visited (x, y) coordinate tuples.
        
        Returns:
            The number of unique self-avoiding walks from the current state.
        """
        # Base case: if no steps are left, we have found one complete walk.
        if steps_left == 0:
            return 1

        count = 0
        x, y = current_pos
        # The four possible moves on a Manhattan (square) lattice.
        moves = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

        for next_pos in moves:
            # If the next position is not already in our path, explore it.
            if next_pos not in visited:
                # Add the new position to the path.
                visited.add(next_pos)
                # Recursively call to continue the walk.
                count += count_walks_recursive(steps_left - 1, next_pos, visited)
                # Backtrack: remove the point to explore other paths from the current position.
                visited.remove(next_pos)
        
        return count

    # Due to the grid's symmetry, we can calculate the number of walks for one
    # starting direction (e.g., right) and multiply the result by 4.
    # The walk starts at (0,0). We fix the first step to be to (1,0).
    # This leaves us needing to count the number of (N-1)-step walks from (1,0).
    start_pos = (1, 0)
    
    # The visited set must initially contain both the origin and the first step's position.
    initial_visited = {(0, 0), start_pos}
    
    # Calculate the number of valid walks for this one starting direction.
    count_one_direction = count_walks_recursive(N - 1, start_pos, initial_visited)

    # The total number of walks is 4 times this count.
    total_walks = 4 * count_one_direction
    
    # Print the final equation with all its numeric components as requested.
    print(f"a({N}) = 4 * {count_one_direction} = {total_walks}")

solve()
