import sys

def solve_saw_count():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a 2D lattice.
    """
    # Set a higher recursion limit for deep recursive calls. n=10 is not too
    # deep, but this is good practice for this type of problem.
    sys.setrecursionlimit(2000)

    memo = {}

    def count_walks(steps_left, x, y, visited):
        """
        Recursively counts the number of self-avoiding walks.
        Uses a tuple of sorted visited points as a key for memoization,
        but symmetry breaking in the main function is the primary optimization.
        """
        # Base case: If no steps are left, we have found one complete walk.
        if steps_left == 0:
            return 1

        count = 0
        # Define the four possible moves on a Manhattan lattice.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]  # Up, Down, Right, Left

        # Explore each possible move from the current position (x, y).
        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy

            # If the next position has not been visited yet...
            if (next_x, next_y) not in visited:
                # ...explore from there.
                visited.add((next_x, next_y))
                count += count_walks(steps_left - 1, next_x, next_y, visited)
                # Backtrack: remove the position to allow exploration of other paths.
                visited.remove((next_x, next_y))
        
        return count

    # The number of steps for the walk.
    n = 10

    # Optimization: Due to grid symmetry, we calculate the walks for one
    # initial direction and multiply by 4.
    # We'll calculate for an initial step to the right: (0,0) -> (1,0).
    # This means we need to find the number of 9-step walks starting from (1,0).
    
    start_pos = (0, 0)
    first_step_pos = (1, 0)
    
    # The path starts with the origin and the first step already visited.
    initial_visited = {start_pos, first_step_pos}
    steps_remaining = n - 1
    
    # Calculate the number of walks for one initial direction.
    walks_one_direction = count_walks(steps_remaining, first_step_pos[0], first_step_pos[1], initial_visited)
    
    # The total number of walks is 4 times this value.
    total_walks = 4 * walks_one_direction
    
    print("Let a(n) be the number of n-step self-avoiding walks on a Manhattan lattice.")
    print(f"To find a({n}), we use a recursive backtracking algorithm with a symmetry optimization.")
    print("Total walks a(n) = 4 * (number of walks starting with a step in one direction, e.g., right).")
    print(f"Number of {n-1}-step walks from (1,0) with (0,0) already visited: {walks_one_direction}")
    print(f"The final equation is: a({n}) = 4 * {walks_one_direction}")
    print(f"Result: a({n}) = {total_walks}")

solve_saw_count()