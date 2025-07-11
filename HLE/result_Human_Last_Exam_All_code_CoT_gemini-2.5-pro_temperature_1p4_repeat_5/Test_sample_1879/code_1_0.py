import sys

# The recursion depth could be up to 4+8=12.
# The default Python recursion limit is usually high enough, but setting it higher is a safe practice for recursive solutions.
sys.setrecursionlimit(2000)

# Define target coordinates and movement constraints
TARGET_X = 4
TARGET_Y = 8
# "four or more" consecutive steps are forbidden, so the maximum allowed is 3.
MAX_CONSECUTIVE_ALLOWED = 3

# Memoization table to store results of subproblems
memo = {}

def count_unique_ways(x, y, last_dir, consecutive):
    """
    Recursively counts the number of unique paths from (x, y) to the target
    with the given constraints.
    """
    # If the path goes out of bounds, it's invalid.
    if x > TARGET_X or y > TARGET_Y:
        return 0

    # If we have reached the target, we found a valid path.
    if x == TARGET_X and y == TARGET_Y:
        return 1

    # Check the memoization table to see if we've already computed this state.
    state = (x, y, last_dir, consecutive)
    if state in memo:
        return memo[state]

    total_ways = 0

    # Option 1: Try to move Right
    if x < TARGET_X:
        if last_dir == 'R':
            # Continuing a sequence of Right moves.
            # We can only do this if we have made fewer than the max allowed consecutive moves.
            if consecutive < MAX_CONSECUTIVE_ALLOWED:
                total_ways += count_unique_ways(x + 1, y, 'R', consecutive + 1)
        else:
            # Starting a new sequence of Right moves.
            total_ways += count_unique_ways(x + 1, y, 'R', 1)

    # Option 2: Try to move Up
    if y < TARGET_Y:
        if last_dir == 'U':
            # Continuing a sequence of Up moves.
            if consecutive < MAX_CONSECUTIVE_ALLOWED:
                total_ways += count_unique_ways(x, y + 1, 'U', consecutive + 1)
        else:
            # Starting a new sequence of Up moves.
            total_ways += count_unique_ways(x, y + 1, 'U', 1)

    # Store the result in the memoization table before returning.
    memo[state] = total_ways
    return total_ways

# A valid path from (0,0) must start with either a Right move or an Up move.
# We calculate the number of paths for each starting move and add them together.

# Calculate paths starting with a Right move. We are now at (1,0).
# The last move was 'R', and the consecutive count is 1.
ways_starting_with_R = count_unique_ways(1, 0, 'R', 1)

# Calculate paths starting with an Up move. We are now at (0,1).
# The last move was 'U', and the consecutive count is 1.
ways_starting_with_U = count_unique_ways(0, 1, 'U', 1)

# The total number of unique ways is the sum of the two cases.
total_unique_ways = ways_starting_with_R + ways_starting_with_U

# Print the final result as an equation
print(f"Number of paths starting with Right: {ways_starting_with_R}")
print(f"Number of paths starting with Up: {ways_starting_with_U}")
print(f"Total number of unique ways is the sum:")
print(f"{ways_starting_with_R} + {ways_starting_with_U} = {total_unique_ways}")
<<<155>>>