import sys

# It's possible to hit the recursion limit for larger T, so we can increase it.
# The chosen limit of 2000 is safe for the search range.
sys.setrecursionlimit(2000)

# Memoization cache to store the results for states (n, p) we've already solved.
# This avoids re-computation and speeds up the process significantly.
# The key is a tuple (n, p), and the value is a boolean (True if winning, False if losing).
memo = {}

def is_winning(n, p):
    """
    Determines if the current state (n tokens remaining, previous move size p)
    is a winning position for the current player using recursion and memoization.
    """
    # Base case: If there are no tokens left, the current player cannot move, so they have lost.
    if n == 0:
        return False
    # Check the cache for a previously computed result to avoid redundant calculations.
    if (n, p) in memo:
        return memo[(n, p)]

    # The current player can take x tokens. The number of tokens taken, x, must be:
    # 1 <= x <= n (must take at least one, and cannot take more than available)
    # x <= 2 * p (the main game rule)
    max_move = min(n, 2 * p)

    # A position is winning if there is at least one move to a losing position.
    # We iterate through all possible moves to find one such move.
    for x in range(1, max_move + 1):
        # The next state for the opponent will be (n - x) tokens, with the previous move being x.
        # If this next state is a losing one for the opponent, then we have found a winning move.
        if not is_winning(n - x, x):
            # We found a winning move. Cache and return True.
            memo[(n, p)] = True
            return True

    # If the loop completes, it means every possible move leads to a winning
    # position for the opponent. Therefore, this current position is a losing one.
    memo[(n, p)] = False
    return False

def find_p2_winning_t(max_t):
    """
    Finds the initial number of tokens T (up to max_t) for which Player 2
    has a winning strategy. This occurs if T is a "losing position" for Player 1.
    """
    p2_winning_numbers = []
    for t in range(2, max_t + 1):
        # Assume T is a losing position for P1 (i.e., P2 wins) until proven otherwise.
        is_p1_losing_pos = True
        
        # Player 1 can make a first move x1, where 1 <= x1 < t.
        # If any of these moves leads to a losing position for Player 2,
        # then Player 1 has a winning strategy, and T is not a solution.
        for x1 in range(1, t):
            # After P1 takes x1 tokens, P2 faces (t - x1) tokens, and the previous move was x1.
            # We check if this is a winning position for P2.
            if not is_winning(t - x1, x1):
                # P1 found a move (x1) that puts P2 in a losing position.
                # Therefore, T is a winning position for P1, so we can stop checking this t.
                is_p1_losing_pos = False
                break
        
        # If after checking all of P1's possible moves, none led to a losing position for P2,
        # then T is indeed a losing position for P1.
        if is_p1_losing_pos:
            p2_winning_numbers.append(t)
            
    return p2_winning_numbers

# Set a limit to search for the numbers. The analysis shows these numbers are Fibonacci numbers.
# A search limit of 35 will find all such numbers up to 34.
search_limit = 35
winning_t = find_p2_winning_t(search_limit)

print("The second player has a winning strategy when the total number of tokens T is one of the following:")
# The problem asks to output each number.
for number in winning_t:
    print(number)
