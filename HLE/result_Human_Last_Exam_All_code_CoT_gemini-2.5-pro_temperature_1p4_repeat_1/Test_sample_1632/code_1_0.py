def count_self_avoiding_walks(n):
    """
    Calculates the number of n-step self-avoiding walks on a 2D square lattice
    using a recursive backtracking algorithm.
    """
    
    # Memoization dictionary (optional, but can speed up for more complex problems)
    # For n=10, direct recursion is fast enough.
    # memo = {}

    def solve(x, y, steps_left, visited):
        """
        Recursively counts valid walks.
        x, y: current coordinates
        steps_left: number of steps remaining
        visited: a set of (x,y) tuples representing the path so far
        """
        # Base case: If no steps are left, we have found one complete valid walk.
        if steps_left == 0:
            return 1
        
        # State for memoization key (if used)
        # key = (x, y, steps_left, tuple(sorted(list(visited))))
        # if key in memo:
        #     return memo[key]

        count = 0
        # Explore the four neighboring positions
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_x, next_y = x + dx, y + dy
            
            # If the next position has not been visited yet
            if (next_x, next_y) not in visited:
                # 1. Add the next position to the path
                visited.add((next_x, next_y))
                
                # 2. Recursively call for the rest of the steps
                count += solve(next_x, next_y, steps_left - 1, visited)
                
                # 3. Backtrack: remove the position to allow other paths to use it
                visited.remove((next_x, next_y))

        # memo[key] = count # Store result if memoization is used
        return count

    # Start the walk from the origin (0,0)
    # The initial path contains just the starting point.
    start_pos = (0, 0)
    initial_visited = {start_pos}
    
    return solve(start_pos[0], start_pos[1], n, initial_visited)

# Define the number of steps
N = 10

# Calculate the result
result = count_self_avoiding_walks(N)

# Print the final equation as requested
print(f"a({N}) = {result}")
