import functools
import sys

# For a walk of length 10, the recursion depth can be significant.
# We increase Python's default recursion limit to prevent potential errors.
sys.setrecursionlimit(2000)

@functools.lru_cache(maxsize=None)
def count_walks(steps_left, current_pos, visited_points):
    """
    Recursively counts the number of n-step self-avoiding walks on a grid.
    This function is memoized using @lru_cache for performance.

    Args:
        steps_left: The number of steps remaining in the walk.
        current_pos: An (x, y) tuple for the current position.
        visited_points: A frozenset of (x, y) tuples of previously visited points.
                        A frozenset is used because it's hashable and can be cached.
    """
    # Base case: If there are no steps left, we have found one complete walk.
    if steps_left == 0:
        return 1

    x, y = current_pos
    total_walks = 0

    # Consider the four possible moves: up, down, left, right.
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    for dx, dy in moves:
        next_pos = (x + dx, y + dy)

        # A walk is self-avoiding, so we can only move to an unvisited square.
        if next_pos not in visited_points:
            # Recursively call the function for the next step.
            # We pass a new frozenset that includes the next position.
            total_walks += count_walks(
                steps_left - 1,
                next_pos,
                visited_points.union({next_pos})
            )

    return total_walks

# The number of steps for the walk.
n = 10

# All walks start at the origin (0, 0).
start_pos = (0, 0)

# The initial set of visited points contains only the starting position.
# We use a frozenset to make it compatible with the cache.
initial_visited = frozenset([start_pos])

# Calculate the total number of n-step self-avoiding walks.
result = count_walks(n, start_pos, initial_visited)

# Print the final result, showing the numbers used in the "equation".
print(f"a({n}) = {result}")