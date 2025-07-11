import sys

# It's possible to hit the recursion limit for larger T, so we increase it.
sys.setrecursionlimit(2000)

# Memoization cache to store the results for states (n, p) we have already computed.
# This drastically speeds up the computation.
memoization_cache = {}

def is_losing_position(n, p):
    """
    Determines if the current game state (n tokens left, p tokens taken previously)
    is a losing position (L-position) for the current player.

    A position is "losing" if every possible move leads to a "winning" position
    for the opponent.
    """
    # If we have computed this state before, return the cached result.
    if (n, p) in memoization_cache:
        return memoization_cache[(n, p)]

    # A player can take x tokens, where 1 <= x <= 2 * p.
    # If they can take all remaining tokens (n), they win immediately.
    # So, if n <= 2 * p, this is a winning position, not a losing one.
    if n <= 2 * p:
        memoization_cache[(n, p)] = False
        return False

    # A position is "winning" if there is at least one move leading to a
    # "losing" position for the opponent.
    # We check all possible moves to see if we can find such a move.
    max_move = min(n, 2 * p)
    for x in range(1, max_move + 1):
        # If the opponent's next position (n - x, x) is a losing one,
        # it means we have found a winning move from the current position.
        if is_losing_position(n - x, x):
            # This means the current position (n, p) is a winning one.
            memoization_cache[(n, p)] = False
            return False

    # If after checking all possible moves, none of them lead to
    # a losing position for the opponent, then this current position
    # must be a losing one.
    memoization_cache[(n, p)] = True
    return True


def find_p2_winning_starts(limit):
    """
    Finds and prints all initial token counts T up to a given limit
    for which the second player (P2) has a winning strategy.

    P2 wins if, for every possible first move by P1, P2 is left in a
    winning position.
    """
    print(f"Searching for T values up to {limit} where P2 wins...")
    p2_win_values = []
    
    for T in range(1, limit + 1):
        p1_can_force_loss = False
        # If T=1, P1 has no valid moves (1 <= x1 < 1 is impossible). P1 loses, so P2 wins.
        if T == 1:
            p2_win_values.append(T)
            continue

        # P1 can take x1 tokens, where 1 <= x1 < T
        for x1 in range(1, T):
            # After P1 takes x1 tokens, the state for P2 is (T - x1, x1).
            # If this state is a losing position for P2, then P1 has found
            # a winning move, and P2 does not have a guaranteed win for this T.
            if is_losing_position(T - x1, x1):
                p1_can_force_loss = True
                break # P1 has a winning move, no need to check other x1
        
        # If after checking all of P1's moves, P1 could not force a loss, then P2 wins.
        if not p1_can_force_loss:
            p2_win_values.append(T)
    
    print("\nThe values of T for which the second player has a winning strategy are:")
    # The theory states these are the Fibonacci numbers. Our code discovers this.
    print(p2_win_values)


# Set a limit for T to check. The calculation can be slow for T > 40.
# A limit of 35 will find the first few numbers in the sequence.
search_limit = 35
find_p2_winning_starts(search_limit)