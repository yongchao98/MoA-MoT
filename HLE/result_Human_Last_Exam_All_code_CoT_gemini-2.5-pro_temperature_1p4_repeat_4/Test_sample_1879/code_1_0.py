import sys

# Increase recursion limit for deep search paths
sys.setrecursionlimit(2000)

# Define the target coordinates and the constraint on consecutive moves.
# The problem forbids 4 or more, so the maximum allowed is 3.
TARGET_R = 4
TARGET_U = 8
MAX_CONSECUTIVE = 3

# Memoization table to store results of subproblems
memo = {}

def count_paths(r, u, last_move, consecutive):
    """
    Recursively counts valid paths from (r, u) to the target.
    
    Args:
        r: Current number of right moves.
        u: Current number of up moves.
        last_move: 'R' or 'U' for the last move made.
        consecutive: Count of consecutive last_move.
        
    Returns:
        The number of unique valid paths from the current state.
    """
    # Use a tuple as a key for the memoization table
    state = (r, u, last_move, consecutive)
    
    # If this state has been computed before, return the stored value
    if state in memo:
        return memo[state]

    # Base case: If we've successfully reached the target, this is one valid path
    if r == TARGET_R and u == TARGET_U:
        return 1
    
    # Start counting ways from this state
    ways = 0
    
    # Option 1: Try to move Right
    # Condition: We have not exceeded the allowed number of right moves (r < 4)
    # AND we are not violating the consecutive move rule.
    if r < TARGET_R and not (last_move == 'R' and consecutive == MAX_CONSECUTIVE):
        # If the last move was 'R', increment consecutive count. Otherwise, reset to 1.
        new_consecutive = consecutive + 1 if last_move == 'R' else 1
        ways += count_paths(r + 1, u, 'R', new_consecutive)

    # Option 2: Try to move Up
    # Condition: We have not exceeded the allowed number of up moves (u < 8)
    # AND we are not violating the consecutive move rule.
    if u < TARGET_U and not (last_move == 'U' and consecutive == MAX_CONSECUTIVE):
        # If the last move was 'U', increment consecutive count. Otherwise, reset to 1.
        new_consecutive = consecutive + 1 if last_move == 'U' else 1
        ways += count_paths(r, u + 1, 'U', new_consecutive)
        
    # Store the computed result in the memo before returning
    memo[state] = ways
    return ways

def solve():
    """
    Calculates the total number of paths and prints the breakdown.
    """
    # Any path must start with either a Right move or an Up move.
    # We can calculate the total number of ways for each starting move.
    
    # Calculate paths starting with a Right move
    ways_starting_R = count_paths(1, 0, 'R', 1)
    
    # Calculate paths starting with an Up move
    ways_starting_U = count_paths(0, 1, 'U', 1)

    # The total is the sum of these two disjoint sets of paths
    total_ways = ways_starting_R + ways_starting_U

    # Print the breakdown as a final equation
    print("The total number of unique ways is the sum of paths starting with 'R' and paths starting with 'U'.")
    print("Final Equation:")
    print(f"{ways_starting_R} (starting with R) + {ways_starting_U} (starting with U) = {total_ways}")
    print("\nTotal unique ways:")
    print(total_ways)

if __name__ == '__main__':
    solve()