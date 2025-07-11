# Target coordinates
TARGET_X = 4
TARGET_Y = 8
# Maximum number of allowed consecutive moves in the same direction
MAX_CONSECUTIVE = 3

# Memoization cache to store results of subproblems
memo = {}

def count_paths(x, y, last_move, consecutive):
    """
    Recursively counts the number of valid paths from (x, y) to the target.

    Args:
        x (int): Current x-coordinate.
        y (int): Current y-coordinate.
        last_move (str): 'R' for right, 'U' for up.
        consecutive (int): Number of consecutive moves in the same direction.

    Returns:
        int: The number of unique valid paths from the current state.
    """
    # Base case: If we have made more than MAX_CONSECUTIVE moves, this path is invalid.
    if consecutive > MAX_CONSECUTIVE:
        return 0

    # Base case: If we have reached the target, we have found one valid path.
    if x == TARGET_X and y == TARGET_Y:
        return 1

    # Base case: If we have overshot the target, this path is invalid.
    if x > TARGET_X or y > TARGET_Y:
        return 0

    # Create a state tuple for memoization
    state = (x, y, last_move, consecutive)
    # If we have already computed this state, return the stored result.
    if state in memo:
        return memo[state]

    # Recursive step: Explore moving right and moving up from the current position.
    
    # Calculate ways if the next move is Right.
    # If the last move was also Right, increment the consecutive counter.
    # Otherwise, reset the counter to 1 for the new direction.
    if last_move == 'R':
        ways = count_paths(x + 1, y, 'R', consecutive + 1)
    else:
        ways = count_paths(x + 1, y, 'R', 1)
        
    # Calculate ways if the next move is Up and add to the total.
    # If the last move was also Up, increment the consecutive counter.
    # Otherwise, reset the counter to 1 for the new direction.
    if last_move == 'U':
        ways += count_paths(x, y + 1, 'U', consecutive + 1)
    else:
        ways += count_paths(x, y + 1, 'U', 1)

    # Store the result in the cache before returning.
    memo[state] = ways
    return ways

if __name__ == "__main__":
    # The starting point is (0,0). The first move can be Right or Up.
    # We can think of this as starting from a "None" state at (0,0).
    # The total number of ways is the sum of paths starting with a Right move
    # and paths starting with an Up move.
    
    ways_starting_R = count_paths(1, 0, 'R', 1)
    ways_starting_U = count_paths(0, 1, 'U', 1)
    
    total_ways = ways_starting_R + ways_starting_U
    
    print("The problem asks for the number of paths from (0,0) to (4,8) with no more than 3 consecutive moves in the same direction.")
    print("This can be calculated by summing the valid paths starting with a 'Right' move and those starting with an 'Up' move.")
    print("\nFinal Equation:")
    print(f"Number of ways = (Ways starting with Right) + (Ways starting with Up)")
    print(f"Number of ways = {ways_starting_R} + {ways_starting_U} = {total_ways}")
    print("\n<<<155>>>")
