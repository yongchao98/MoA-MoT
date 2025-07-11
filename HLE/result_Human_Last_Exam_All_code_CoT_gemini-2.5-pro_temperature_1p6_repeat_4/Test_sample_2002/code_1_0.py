import sys

# Increase the recursion limit for checking larger values of T, if necessary.
# The default limit is usually 1000, which should be sufficient for T < 100.
sys.setrecursionlimit(2000)

# Memoization cache to store the results for (r, l) pairs.
# This is crucial for performance, as many subproblems are repeated.
memo = {}

def is_losing_pos(r, l):
    """
    Determines if the current game position is a losing one for the current player.
    A position is defined by (r) remaining tokens and the (l) tokens taken in the last move.
    
    A position is losing if all possible moves lead to a winning position for the opponent.
    Returns:
        True if the position (r, l) is a losing one, False otherwise.
    """
    # Base case: If there are no tokens left, you cannot move, so you have lost.
    if r == 0:
        return True
    
    # Check if we have already computed the result for this state.
    if (r, l) in memo:
        return memo[(r, l)]

    # The maximum number of tokens a player can take is 2 times the last move,
    # but not more than the remaining tokens.
    max_move = min(r, 2 * l)
    
    # Iterate through all possible moves 'x'.
    for x in range(1, max_move + 1):
        # The next state for the opponent will be (r - x, x).
        # If we can find a move 'x' that puts the opponent in a losing position,
        # it means the current position is a WINNING one.
        if is_losing_pos(r - x, x):
            # We found a winning move. So, this position is NOT losing.
            memo[(r, l)] = False
            return False

    # If the loop completes, it means no move 'x' leads to a losing position for the opponent.
    # Therefore, all our moves lead to a winning position for them.
    # This means the current position (r, l) is a LOSING one.
    memo[(r, l)] = True
    return True

def is_p2_win(t):
    """
    Determines if Player 2 has a winning strategy for an initial pile of T tokens.
    
    Player 2 wins if for EVERY possible first move by Player 1, Player 2 can force a win.
    """
    # Player 1's first move 'n1' can be any number of tokens from 1 to t-1.
    for n1 in range(1, t):
        # After P1 takes n1 tokens, r = t - n1 tokens are left.
        # It becomes P2's turn, and the last move was n1.
        # P2 is now at position (t - n1, n1).
        # If this position is a losing one for P2, then P1 has a winning strategy.
        if is_losing_pos(t - n1, n1):
            # P1 found a move (n1) that puts P2 in a losing position.
            # This means T is a win for P1, not P2.
            return False

    # If the loop completes, it means P1 cannot force P2 into a losing position.
    # Therefore, P2 has a winning strategy for this T.
    return True

def find_p2_winning_numbers(limit):
    """
    Finds and prints all T values up to a given limit for which Player 2 wins.
    """
    print(f"Searching for T values up to {limit} where the second player wins...")
    p2_win_numbers = []
    # The game is non-trivial for T > 1. T=1 is a degenerate case where P1 has no moves.
    # Let's check from T=1.
    for T in range(1, limit + 1):
        if is_p2_win(T):
            p2_win_numbers.append(T)
            
    print("The second player has a winning strategy for the following T values:")
    print(p2_win_numbers)
    
    fibs = [1, 2]
    while fibs[-1] <= limit:
        fibs.append(fibs[-1] + fibs[-2])
    print("\nFor comparison, the Fibonacci numbers are:")
    # Print Fibonacci numbers up to the limit
    print([f for f in fibs if f <= limit])

if __name__ == '__main__':
    # We check for T up to 90. F_11=89.
    find_p2_winning_numbers(90)