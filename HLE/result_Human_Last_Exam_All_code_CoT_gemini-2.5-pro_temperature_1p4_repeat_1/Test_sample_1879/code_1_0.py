import sys

# It's good practice to increase the recursion limit for recursive path-counting problems,
# although the default may be sufficient for this specific grid size.
sys.setrecursionlimit(2000)

# The destination coordinates
TARGET_X = 4
TARGET_Y = 8

# A cache to store the results of subproblems (memoization)
memo = {}

def count_unique_ways(x, y, last_move, consecutive_steps):
    """
    Calculates the number of unique paths from (x, y) to the target,
    respecting the constraint on consecutive moves.

    Args:
        x: The current x-coordinate.
        y: The current y-coordinate.
        last_move: The direction of the last move ('R' for Right, 'U' for Up, or 'S' for Start).
        consecutive_steps: The number of consecutive steps taken in the last_move direction.

    Returns:
        The number of unique ways from the current state to the target.
    """
    # Base case: If we have reached the target, we have found one valid path.
    if x == TARGET_X and y == TARGET_Y:
        return 1

    # Base case: If we have moved past the target, this path is invalid.
    if x > TARGET_X or y > TARGET_Y:
        return 0

    # Create a unique key for the current state to use with the memoization cache.
    state = (x, y, last_move, consecutive_steps)
    if state in memo:
        return memo[state]

    total_ways = 0

    # --- Option 1: Try moving RIGHT ---
    # We can move right as long as the last move wasn't a run of 3 'R's.
    if last_move == 'R':
        # If the last move was 'R', we can continue the run if it's not too long.
        if consecutive_steps < 3:
            total_ways += count_unique_ways(x + 1, y, 'R', consecutive_steps + 1)
    else:
        # If the last move was 'U' or 'S' (start), we can always start a new 'R' run.
        total_ways += count_unique_ways(x + 1, y, 'R', 1)

    # --- Option 2: Try moving UP ---
    # We can move up as long as the last move wasn't a run of 3 'U's.
    if last_move == 'U':
        # If the last move was 'U', we can continue the run if it's not too long.
        if consecutive_steps < 3:
            total_ways += count_unique_ways(x, y + 1, 'U', consecutive_steps + 1)
    else:
        # If the last move was 'R' or 'S' (start), we can always start a new 'U' run.
        total_ways += count_unique_ways(x, y + 1, 'U', 1)

    # Cache the result for the current state before returning it.
    memo[state] = total_ways
    return total_ways

# Main execution:
# Start the recursion from (0,0) with an initial 'Start' state.
total_paths = count_unique_ways(0, 0, 'S', 0)

# The total number of paths is the sum of paths that begin with 'R' and those that begin with 'U'.
# Our function calculates these subproblems, and the results are stored in the memo cache.
ways_starting_R = memo.get((1, 0, 'R', 1), 0)
ways_starting_U = memo.get((0, 1, 'U', 1), 0)

print(f"To get from (0,0) to ({TARGET_X},{TARGET_Y}) with at most 3 consecutive steps in one direction:")
print(f"Number of paths starting with a step to the Right: {ways_starting_R}")
print(f"Number of paths starting with a step Up: {ways_starting_U}")
print(f"The total number of unique ways is the sum of these two possibilities.")
print(f"Final Equation: {ways_starting_R} + {ways_starting_U} = {total_paths}")