import collections

# Target coordinates
TARGET_R = 4
TARGET_U = 8

# Maximum number of consecutive steps in the same direction.
# The constraint is "cannot move four or more", so the max allowed is 3.
MAX_CONSECUTIVE = 3

# Memoization table to store results of computed states
memo = {}

def solve(r, u, last_move, consecutive):
    """
    Calculates the number of unique paths from (r, u) to the target.

    Args:
        r: current right steps (x-coordinate)
        u: current up steps (y-coordinate)
        last_move: 0 for Right, 1 for Up
        consecutive: number of consecutive steps for last_move
    Returns:
        The number of unique paths from the current state.
    """
    # If the state is already computed, return the stored value
    state = (r, u, last_move, consecutive)
    if state in memo:
        return memo[state]

    # Base case: If we have reached the destination, we found a valid path
    if r == TARGET_R and u == TARGET_U:
        return 1

    # Base case: If we have overshot the destination, this path is invalid
    if r > TARGET_R or u > TARGET_U:
        return 0

    total_ways = 0

    # --- Try to move Right ---
    if r < TARGET_R:
        # If the last move was also Right
        if last_move == 0:
            # We can continue if we haven't reached the max consecutive steps
            if consecutive < MAX_CONSECUTIVE:
                total_ways += solve(r + 1, u, 0, consecutive + 1)
        # If the last move was Up, we can switch to Right
        else:  # last_move == 1
            total_ways += solve(r + 1, u, 0, 1)

    # --- Try to move Up ---
    if u < TARGET_U:
        # If the last move was also Up
        if last_move == 1:
            # We can continue if we haven't reached the max consecutive steps
            if consecutive < MAX_CONSECUTIVE:
                total_ways += solve(r, u + 1, 1, consecutive + 1)
        # If the last move was Right, we can switch to Up
        else:  # last_move == 0
            total_ways += solve(r, u + 1, 1, 1)

    # Store the result in the memoization table and return it
    memo[state] = total_ways
    return total_ways

def main():
    """
    Main function to solve the problem.
    The path starts at (0,0). The first step can be either Right or Up.
    We calculate the number of paths for each starting move and sum them up.
    """
    # Calculate paths starting with a Right move from (0,0) to (1,0)
    # State at (1,0) is: r=1, u=0, last_move=Right(0), consecutive=1
    ways_starting_with_R = solve(1, 0, 0, 1)

    # Calculate paths starting with an Up move from (0,0) to (0,1)
    # State at (0,1) is: r=0, u=1, last_move=Up(1), consecutive=1
    ways_starting_with_U = solve(0, 1, 1, 1)

    # The total number of ways is the sum of the two possibilities
    total_ways = ways_starting_with_R + ways_starting_with_U

    print("To find the total number of unique ways, we can sum the ways for paths that start with a 'Right' move and those that start with an 'Up' move.")
    print(f"Number of ways starting with 'Right': {ways_starting_with_R}")
    print(f"Number of ways starting with 'Up': {ways_starting_with_U}")
    print("\nThe final equation is:")
    print(f"{ways_starting_with_R} + {ways_starting_with_U} = {total_ways}")

if __name__ == "__main__":
    main()