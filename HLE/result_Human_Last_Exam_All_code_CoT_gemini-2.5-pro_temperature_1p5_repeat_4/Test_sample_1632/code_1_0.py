import sys

def solve_saw_count():
    """
    Calculates the number of n-step self-avoiding walks on a Manhattan lattice.
    """
    # The length of the walk
    n = 10
    
    # Increase recursion limit for deep searches, although n=10 is manageable.
    sys.setrecursionlimit(n * 4 + 10)

    # memo = {} # Memoization is tricky due to the 'visited' set state.
    # We will use a direct recursive count.

    def count_walks(x, y, steps_left, visited):
        """
        Recursively counts the number of self-avoiding walks.
        x, y: current coordinates
        steps_left: number of steps remaining
        visited: a set of (x,y) tuples of visited points
        """
        # Base case: if no steps are left, we have found one valid walk.
        if steps_left == 0:
            return 1

        total = 0
        
        # Explore neighbors: Right(1,0), Left(-1,0), Up(0,1), Down(0,-1)
        moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]

        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy
            if (next_x, next_y) not in visited:
                # Add the new point to the path
                visited.add((next_x, next_y))
                # Recurse from the new point
                total += count_walks(next_x, next_y, steps_left - 1, visited)
                # Backtrack: remove the point to explore other paths
                visited.remove((next_x, next_y))
        
        return total

    # --- Main Calculation ---
    print(f"Calculating a({n}), the number of {n}-step self-avoiding walks.")
    
    if n == 0:
        result = 1
        print("a(0) = 1")
    else:
        # Optimization: Use 4-fold symmetry of the lattice.
        # We calculate the number of walks for one initial direction (e.g., Right)
        # and multiply the result by 4.
        
        # This means the walk starts at (0,0), and the first step is to (1,0).
        # We start our recursive search from (1,0) with n-1 steps left.
        start_x, start_y = 1, 0
        steps_left = n - 1
        # The initial visited set includes the origin and the first step.
        initial_visited = {(0, 0), (1, 0)}

        walks_one_direction = count_walks(start_x, start_y, steps_left, initial_visited)
        total_walks = 4 * walks_one_direction

        # Output the final equation as requested
        print(f"Due to symmetry, we can find the number of walks starting in one direction and multiply by 4.")
        print(f"The number of walks starting with a step to the right is: {walks_one_direction}")
        print(f"a({n}) = 4 * {walks_one_direction}")
        print(f"a({n}) = {total_walks}")

solve_saw_count()
<<<259140>>>