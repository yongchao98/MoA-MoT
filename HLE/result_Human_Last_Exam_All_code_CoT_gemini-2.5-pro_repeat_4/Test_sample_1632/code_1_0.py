import collections

def solve():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a 2D square lattice.
    """
    n = 10
    # Memoization cache to store results of subproblems.
    # The key is a tuple: (current_x, current_y, steps_left, frozenset_of_visited_points)
    memo = {}

    def count_walks(x, y, steps_left, visited):
        """
        Recursively counts the number of self-avoiding walks.
        """
        # Base case: If no steps are left, we have found one valid walk.
        if steps_left == 0:
            return 1

        # Create a key for the memoization cache.
        # The set of visited points must be converted to an immutable frozenset to be used as a dictionary key.
        state = (x, y, steps_left, frozenset(visited))
        if state in memo:
            return memo[state]

        total = 0
        # Explore the four possible directions: up, down, left, right.
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nx, ny = x + dx, y + dy

            # If the next point has not been visited, explore from there.
            if (nx, ny) not in visited:
                new_visited = visited.copy()
                new_visited.add((nx, ny))
                total += count_walks(nx, ny, steps_left - 1, new_visited)

        # Cache the result for the current state and return it.
        memo[state] = total
        return total

    # Initial setup:
    # Start at the origin (0,0).
    # Use symmetry: calculate for one initial direction and multiply by 4.
    # We take the first step from (0,0) to (1,0).
    # This leaves n-1 = 9 steps to be calculated.
    start_x, start_y = 0, 0
    first_step_x, first_step_y = 1, 0
    
    initial_visited = {(start_x, start_y), (first_step_x, first_step_y)}
    steps_remaining = n - 1

    # Calculate the number of walks for one initial direction.
    walks_one_direction = count_walks(first_step_x, first_step_y, steps_remaining, initial_visited)

    # The total number of walks is 4 times this result.
    total_walks = 4 * walks_one_direction

    print(f"For a walk of length n={n}:")
    print(f"Number of paths starting with a single direction (e.g., right): {walks_one_direction}")
    print(f"Total number of self-avoiding walks a({n}) = 4 * {walks_one_direction} = {total_walks}")

solve()
<<<259140>>>