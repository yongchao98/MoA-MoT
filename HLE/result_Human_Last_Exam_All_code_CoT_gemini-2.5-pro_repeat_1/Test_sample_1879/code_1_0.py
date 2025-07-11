import sys

# Increase the recursion limit for potentially deep recursive calls.
sys.setrecursionlimit(2000)

# Use a dictionary for memoization to cache the results of subproblems.
memo = {}
TARGET_R, TARGET_U = 4, 8

def solve(r, u, is_last_r, consecutive):
    """
    Calculates the number of valid paths from (r, u) to the target.

    Args:
      r: The number of Right moves made so far.
      u: The number of Up moves made so far.
      is_last_r: A boolean that is True if the last move was Right, False if Up.
      consecutive: The number of consecutive moves made in the last direction.
    """
    # Define the state as a tuple to use as a key in the memoization cache.
    state = (r, u, is_last_r, consecutive)
    if state in memo:
        return memo[state]

    # Base case: If we have reached the destination (4, 8), we found one valid path.
    if r == TARGET_R and u == TARGET_U:
        return 1

    # Base case: If we have moved beyond the destination, this path is invalid.
    if r > TARGET_R or u > TARGET_U:
        return 0

    count = 0

    # Option 1: Try to move Right.
    # This is possible if the number of right moves is less than 4.
    if r < TARGET_R:
        if is_last_r:
            # If the last move was also Right, we can only continue if the
            # consecutive count is less than 3.
            if consecutive < 3:
                count += solve(r + 1, u, True, consecutive + 1)
        else:
            # If the last move was Up, we can start a new sequence of Right moves.
            count += solve(r + 1, u, True, 1)

    # Option 2: Try to move Up.
    # This is possible if the number of up moves is less than 8.
    if u < TARGET_U:
        if not is_last_r:
            # If the last move was also Up, we can only continue if the
            # consecutive count is less than 3.
            if consecutive < 3:
                count += solve(r, u + 1, False, consecutive + 1)
        else:
            # If the last move was Right, we can start a new sequence of Up moves.
            count += solve(r, u + 1, False, 1)

    # Cache the result for the current state before returning it.
    memo[state] = count
    return count

# Any path must start with either a Right or an Up move.
# Calculate paths starting with a Right move (from (0,0) to (1,0)).
ways_starting_with_r = solve(1, 0, True, 1)

# Calculate paths starting with an Up move (from (0,0) to (0,1)).
ways_starting_with_u = solve(0, 1, False, 1)

# The total number of unique paths is the sum of these two disjoint sets of paths.
total_ways = ways_starting_with_r + ways_starting_with_u

# Print the final result in the requested equation format.
print(f"{ways_starting_with_r} + {ways_starting_with_u} = {total_ways}")
<<<254>>>