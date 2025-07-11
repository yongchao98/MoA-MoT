import sys

# Increase the recursion limit for deep recursive calls.
# The default limit might be too low for the 8 'Up' steps.
sys.setrecursionlimit(2000)

# Target coordinates
TARGET_R = 4
TARGET_U = 8
# Constraint: Cannot move four or more consecutive steps in the same direction.
# This means the maximum allowed consecutive steps is 3.
MAX_CONSECUTIVE = 3

# Memoization table to store results of subproblems
memo = {}

def count_ways(r, u, last_move, consecutive_steps):
    """
    Recursively counts the number of unique paths from (r, u) to the target,
    respecting the movement constraints.
    """
    # Base case: If we have reached the destination, we found a valid path.
    if r == TARGET_R and u == TARGET_U:
        return 1
    
    # Base case: If we have moved past the destination, this path is invalid.
    if r > TARGET_R or u > TARGET_U:
        return 0

    # Create a unique key for the current state for memoization
    state = (r, u, last_move, consecutive_steps)
    if state in memo:
        return memo[state]

    ways = 0
    # --- Try to move RIGHT ---
    if r < TARGET_R:
        if last_move == 'R':
            # Continue moving right if we haven't exceeded the consecutive limit
            if consecutive_steps < MAX_CONSECUTIVE:
                ways += count_ways(r + 1, u, 'R', consecutive_steps + 1)
        else:
            # If the last move was 'U', we can start a new sequence of 'R' moves
            ways += count_ways(r + 1, u, 'R', 1)
            
    # --- Try to move UP ---
    if u < TARGET_U:
        if last_move == 'U':
            # Continue moving up if we haven't exceeded the consecutive limit
            if consecutive_steps < MAX_CONSECUTIVE:
                ways += count_ways(r, u + 1, 'U', consecutive_steps + 1)
        else:
            # If the last move was 'R', we can start a new sequence of 'U' moves
            ways += count_ways(r, u + 1, 'U', 1)

    # Store the result in the memoization table and return it
    memo[state] = ways
    return ways

# The first step from (0,0) can be either Right or Up.
# We calculate the number of paths for each starting move.

# Case 1: First move is Right. We start the recursion from (1, 0).
ways_starting_with_R = count_ways(1, 0, 'R', 1)

# Case 2: First move is Up. We start the recursion from (0, 1).
# The memoization table already contains results from the first call that can be reused.
ways_starting_with_U = count_ways(0, 1, 'U', 1)

# The total number of ways is the sum of the two cases.
total_ways = ways_starting_with_R + ways_starting_with_U

# Print the final result as an equation
print(f"Number of paths starting with Right: {ways_starting_with_R}")
print(f"Number of paths starting with Up: {ways_starting_with_U}")
print(f"{ways_starting_with_R} + {ways_starting_with_U} = {total_ways}")
