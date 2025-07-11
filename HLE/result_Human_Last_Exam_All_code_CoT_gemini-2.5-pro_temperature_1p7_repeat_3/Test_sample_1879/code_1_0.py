import sys

# It's good practice to increase the recursion limit for potentially deep recursion,
# although for a 4x8 grid, it's not strictly necessary.
sys.setrecursionlimit(2000)

# Target coordinates define the number of Right and Up moves needed.
R_TARGET, U_TARGET = 4, 8

# Memoization table to store results of computed states.
memo = {}

def solve(r, u, is_last_move_right, consecutive_count):
    """
    Calculates the number of valid paths from state (r, u) to the target.
    A state is defined by the number of right/up moves, the last move direction,
    and the count of consecutive moves in that direction.
    """
    
    # Base case: If we've reached the destination, we found one valid path.
    if r == R_TARGET and u == U_TARGET:
        return 1
    
    # Base case: If we've moved past the destination, this path is invalid.
    if r > R_TARGET or u > U_TARGET:
        return 0

    # Create a unique key for the current state to use with the memoization table.
    state = (r, u, is_last_move_right, consecutive_count)
    
    # If this state has been computed before, return the stored result.
    if state in memo:
        return memo[state]

    ways = 0
    
    # Option 1: Try to move Right. This is possible only if r < R_TARGET.
    if r < R_TARGET:
        # Case A: The last move was also Right.
        if is_last_move_right:
            # We can only move Right again if the consecutive count is less than 3.
            if consecutive_count < 3:
                ways += solve(r + 1, u, True, consecutive_count + 1)
        # Case B: The last move was Up (or this is the first move).
        else:
            # We can always switch direction. The new consecutive count becomes 1.
            ways += solve(r + 1, u, True, 1)

    # Option 2: Try to move Up. This is possible only if u < U_TARGET.
    if u < U_TARGET:
        # Case A: The last move was also Up.
        if not is_last_move_right:
            # We can only move Up again if the consecutive count is less than 3.
            if consecutive_count < 3:
                ways += solve(r, u + 1, False, consecutive_count + 1)
        # Case B: The last move was Right.
        else:
            # We can always switch direction. The new consecutive count becomes 1.
            ways += solve(r, u + 1, False, 1)

    # Store the result for the current state in the memoization table and return it.
    memo[state] = ways
    return ways

# ---- Main Execution ----

# Since we start at (0,0), there is no prior move. We calculate the total
# paths by summing the possibilities for the first move.

# Calculate paths that start with a move to the Right.
# The state is at (1,0), last move was Right, with a consecutive count of 1.
ways_starting_with_R = solve(1, 0, True, 1)

# Calculate paths that start with a move Up.
# The state is at (0,1), last move was Up, with a consecutive count of 1.
ways_starting_with_U = solve(0, 1, False, 1)

# The total number of unique ways is the sum of these two scenarios.
total_ways = ways_starting_with_R + ways_starting_with_U

# Print the breakdown of the final calculation as requested.
print(f"Number of ways starting with a Right move: {ways_starting_with_R}")
print(f"Number of ways starting with an Up move: {ways_starting_with_U}")
print(f"Total unique ways = {ways_starting_with_R} + {ways_starting_with_U} = {total_ways}")
<<<71>>>