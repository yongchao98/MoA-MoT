import sys

# It's a recursive problem, so we might hit the recursion limit.
# Set a higher limit, though for MAX_T=100 it's not strictly necessary.
sys.setrecursionlimit(2000)

# A dictionary to store the results of (n, k) states to avoid re-computation.
memo = {}

def is_L(n, k):
    """
    Determines if the game state (n, k) is a losing position for the current player.
    n: number of tokens remaining.
    k: number of tokens taken by the previous player.
    Returns True if the position is losing, False otherwise.
    """
    if n == 0:
        # If there are no tokens left, the current player has no moves and has lost.
        # This is by definition a losing position.
        return True
    if (n, k) in memo:
        return memo[(n, k)]

    # The maximum number of tokens the current player can take.
    # The move must be positive, so x >= 1.
    limit = min(n, 2 * k)

    # A position is a Winning-position if there exists at least one move
    # that leads to a Losing-position for the opponent.
    for x in range(1, limit + 1):
        # The opponent will be in state (n - x, x).
        # If is_L(n - x, x) is True, it means we can force the opponent into a losing position.
        if is_L(n - x, x):
            # This means the current position (n, k) is a Winning-position.
            memo[(n, k)] = False
            return False

    # If no move leads to a losing position for the opponent, it means all moves
    # lead to a winning position for them.
    # Therefore, the current position (n, k) is a Losing-position.
    memo[(n, k)] = True
    return True

def find_p2_winning_T(max_T):
    """
    Finds all values of T up to max_T for which the second player has a winning strategy.
    """
    p2_win_positions = []
    for T in range(1, max_T + 1):
        # Assume T is a P2-win position until we find a counter-example.
        is_p2_win = True

        # P1 can take x1 tokens, where 1 <= x1 < T.
        # If P1 can make any move that puts P2 in a losing position, T is a P1-win.
        for x1 in range(1, T):
            # P2 is in state (T - x1, x1).
            # If this is a losing position for P2, P1 has found a winning move.
            if is_L(T - x1, x1):
                is_p2_win = False
                break
        
        # If after checking all of P1's moves, none could force P2 into a
        # losing position, then T is a P2-win position.
        # The edge case T=1 is handled correctly: range(1,1) is empty, so is_p2_win remains True.
        if is_p2_win:
            p2_win_positions.append(T)
            
    return p2_win_positions

if __name__ == "__main__":
    MAX_T = 100
    print(f"Searching for T values up to {MAX_T} where the second player has a winning strategy...")
    
    winning_Ts = find_p2_winning_T(MAX_T)
    
    print("\nThe second player has a winning strategy for the following initial token counts T:")
    # We print the numbers in the list as requested.
    for t in winning_Ts:
        print(t, end=' ')
    print("\n")
    print("This sequence corresponds to the Fibonacci numbers (F_n where n >= 2, using the convention F_1=1, F_2=1).")
    
    fib = [1, 2]
    while fib[-1] < MAX_T:
        fib.append(fib[-1] + fib[-2])
    # Remove last if it exceeds MAX_T
    if fib[-1] > MAX_T:
        fib.pop()
    print("For comparison, the Fibonacci numbers are:", fib)
