import sys

# The grid dimensions and target coordinates are small enough that
# increasing the recursion limit is not strictly necessary, but it's good practice
# for recursive solutions to avoid potential RecursionError.
sys.setrecursionlimit(2000)

# Memoization table (dictionary) to store results of solved subproblems.
memo = {}
# The destination coordinates.
target_x, target_y = 4, 8

def count_ways(x, y, last_dir, consecutive_len):
    """
    Calculates the number of unique valid paths from coordinates (x, y)
    to the target (4, 8) using recursion with memoization.

    The state is defined by:
    x (int): The current x-coordinate (number of Right steps taken).
    y (int): The current y-coordinate (number of Up steps taken).
    last_dir (int): The direction of the last move. We use 0 for Right and 1 for Up.
                    A value of -1 indicates the starting point with no prior moves.
    consecutive_len (int): The number of consecutive moves made in last_dir.
    """
    # Base Case 1: If we have reached the destination, we've found one valid path.
    if x == target_x and y == target_y:
        return 1

    # Base Case 2: If we have moved past the destination, this path is invalid.
    if x > target_x or y > target_y:
        return 0

    # Create a unique tuple to represent the current state for memoization.
    state = (x, y, last_dir, consecutive_len)
    # If this state has been computed before, return the stored result.
    if state in memo:
        return memo[state]

    # This variable will accumulate the number of ways from the current state.
    total_ways = 0

    # --- Recursive Step ---

    # Option 1: Try to move Right.
    # A move is forbidden if it's the 4th consecutive move in the same direction.
    # So, we can move right if the last move wasn't 'Right' for 3 consecutive steps.
    if not (last_dir == 0 and consecutive_len == 3):
        # If the last move was also Right, we increment the consecutive length.
        # Otherwise, a new sequence of 'Right' moves starts, so the length is 1.
        new_len = consecutive_len + 1 if last_dir == 0 else 1
        total_ways += count_ways(x + 1, y, 0, new_len)

    # Option 2: Try to move Up.
    # Similarly, we can move Up if we haven't already made 3 consecutive Up moves.
    if not (last_dir == 1 and consecutive_len == 3):
        # If the last move was also Up, increment the length. Otherwise, reset to 1.
        new_len = consecutive_len + 1 if last_dir == 1 else 1
        total_ways += count_ways(x, y + 1, 1, new_len)

    # Store the calculated result in our memoization table before returning.
    memo[state] = total_ways
    return total_ways

# We start at point A(0,0). We need to take 4 steps Right and 8 steps Up.
# The constraint is that no more than 3 consecutive steps can be in the same direction.
# The total number of ways can be represented by the equation:
# Total = W(0,0), where W(x,y) is the number of valid paths from (x,y) to (4,8).
# This is solved by the recursive function call below.

# Initial call starts at (0, 0), with no previous move (last_dir=-1, consecutive_len=0).
final_result = count_ways(0, 0, -1, 0)

print(f"The number of unique ways to walk from point A(0,0) to B(4,8) is calculated based on a recursive formula.")
print(f"Let W(x, y, dir, len) be the number of ways from (x, y).")
print(f"W(x,y,...) = W_Right + W_Up, where each term is only added if the move is valid.")
print(f"A move is invalid if it is the 4th consecutive move in the same direction.")
print(f"The starting point is A(0,0) and the destination is B(4,8).")
print(f"Total ways = {final_result}")

<<<155>>>