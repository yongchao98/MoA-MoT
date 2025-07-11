import collections

# Target coordinates
TARGET_X = 4
TARGET_Y = 8

# The constraint is "four or more", which means a maximum of 3 consecutive are allowed.
MAX_CONSECUTIVE = 3

# Memoization table to store results of subproblems
memo = {}

def count_paths(x, y, last_dir, consecutive):
    """
    Recursively counts the number of valid paths from (x, y) to the target.
    
    Args:
        x (int): current x-coordinate.
        y (int): current y-coordinate.
        last_dir (int): Direction of the last move (0 for Right, 1 for Up).
        consecutive (int): Number of consecutive steps in the last_dir.
        
    Returns:
        int: The number of unique valid paths from the current state to the target.
    """
    # Base Case: Successfully reached the target
    if x == TARGET_X and y == TARGET_Y:
        return 1
        
    # Base Case: Moved out of bounds
    if x > TARGET_X or y > TARGET_Y:
        return 0
        
    # Memoization: Check if this state has been computed before
    state = (x, y, last_dir, consecutive)
    if state in memo:
        return memo[state]
        
    # Recursive Step: Calculate ways by moving Right or Up
    total_ways = 0
    
    # Option 1: Move Right
    if x < TARGET_X:
        # If the last move was also Right, we can only continue if we are below the consecutive limit
        if last_dir == 0 and consecutive < MAX_CONSECUTIVE:
            total_ways += count_paths(x + 1, y, 0, consecutive + 1)
        # If the last move was Up, we can always start a new sequence of Right moves
        elif last_dir == 1:
            total_ways += count_paths(x + 1, y, 0, 1)
            
    # Option 2: Move Up
    if y < TARGET_Y:
        # If the last move was also Up, we can only continue if we are below the consecutive limit
        if last_dir == 1 and consecutive < MAX_CONSECUTIVE:
            total_ways += count_paths(x, y + 1, 1, consecutive + 1)
        # If the last move was Right, we can always start a new sequence of Up moves
        elif last_dir == 0:
            total_ways += count_paths(x, y + 1, 1, 1)

    # Store the result in the memoization table before returning
    memo[state] = total_ways
    return total_ways

# The total number of paths is the sum of paths starting with a Right move
# and paths starting with an Up move.

# Calculate paths starting with a Right move (from (0,0) to (1,0))
# The state at (1,0) is: arrived from Right (0), in a sequence of 1.
ways_starting_with_right = count_paths(1, 0, 0, 1)

# Clear memoization table for the next independent calculation.
memo.clear() 

# Calculate paths starting with an Up move (from (0,0) to (0,1))
# The state at (0,1) is: arrived from Up (1), in a sequence of 1.
ways_starting_with_up = count_paths(0, 1, 1, 1)

# The total is the sum of these two disjoint sets of paths.
total_unique_ways = ways_starting_with_right + ways_starting_with_up

print("This problem can be broken down into two main scenarios based on the first step:")
print(f"1. Number of ways starting with a Right move: {ways_starting_with_right}")
print(f"2. Number of ways starting with an Up move: {ways_starting_with_up}")
print("\nThe final equation for the total number of unique ways is the sum of these two scenarios:")
print(f"{ways_starting_with_right} + {ways_starting_with_up} = {total_unique_ways}")
print(f"\nTotal unique ways from A(0,0) to B(4,8) are: {total_unique_ways}")

<<<119>>>