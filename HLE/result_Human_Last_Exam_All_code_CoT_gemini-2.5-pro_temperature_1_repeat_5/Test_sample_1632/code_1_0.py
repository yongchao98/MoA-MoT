import sys

def solve_saw_count():
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice.
    This script is specifically written to compute a(10).
    """
    # Set the desired number of steps for the walk.
    n_steps = 10
    
    # Increase the recursion limit for safety, although n=10 should not exceed the default.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), n_steps + 5))

    # A dictionary for memoization could be used but is complex due to the `visited` set.
    # For n=10, a direct recursive approach is feasible.
    memo = {}

    def count_walks(steps, x, y, visited):
        """
        Recursively counts the number of self-avoiding walks.

        Args:
            steps (int): The number of steps remaining in the walk.
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            visited (frozenset): A frozenset of (x, y) tuples representing visited points.
                                 Using frozenset allows for memoization.
        """
        # Base case: If there are no steps left, we have successfully found one valid path.
        if steps == 0:
            return 1

        # Memoization: Check if we have already computed the result for this state.
        state = (steps, x, y, visited)
        if state in memo:
            return memo[state]

        count = 0
        # Define the four possible moves on a Manhattan lattice: Up, Down, Left, Right.
        moves = [(0, 1), (0, -1), (-1, 0), (1, 0)]

        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy

            # Check if the next point has already been visited.
            if (next_x, next_y) not in visited:
                # If the point is new, continue the walk from there.
                # A new visited set is created for the next recursive call.
                new_visited = visited.union({(next_x, next_y)})
                count += count_walks(steps - 1, next_x, next_y, new_visited)

        # Store the result in the memoization table and return it.
        memo[state] = count
        return count

    # Handle the trivial case of a 0-step walk.
    if n_steps == 0:
        result = 1
    else:
        # Optimization: Utilize the symmetry of the lattice.
        # The number of walks starting with a step in any of the 4 directions is the same.
        # We calculate the number of walks for one initial move (e.g., to (0,1))
        # and multiply by 4. This first step reduces n_steps by 1.
        
        # The walk starts at (0,0). The first step is to (0,1).
        # So, the visited set for the recursive call starts with both points.
        initial_visited = frozenset([(0, 0), (0, 1)])
        
        # Calculate walks for one direction (n-1 steps remaining from (0,1)).
        walks_one_direction = count_walks(n_steps - 1, 0, 1, initial_visited)
        
        # Total number of walks is 4 times the count for one direction.
        result = 4 * walks_one_direction

    # Print the final result, showing each part of the equation as requested.
    print("a", "(", n_steps, ")", "=", result)


solve_saw_count()