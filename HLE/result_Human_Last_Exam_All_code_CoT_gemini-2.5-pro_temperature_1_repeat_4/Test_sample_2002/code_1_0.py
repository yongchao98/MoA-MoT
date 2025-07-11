# A dictionary to store the results of subproblems to avoid recomputation.
memo = {}

def is_winning_state(n, m):
    """
    Recursively determines if a state (n tokens, max_move m) is winning for the current player.
    
    A position is winning if there exists a move to a losing position for the opponent.
    A position is losing if all moves lead to a winning position for the opponent.
    """
    # Base case: If the player can take all remaining tokens, they win.
    if n <= m:
        return True
        
    # Check if we have already computed the result for this state.
    if (n, m) in memo:
        return memo[(n, m)]

    # Iterate through all possible moves (take y tokens).
    for y in range(1, m + 1):
        # The next state for the opponent will be (n-y tokens, max_move 2*y).
        # If this next state is a losing one for the opponent, then the current state is winning.
        if not is_winning_state(n - y, 2 * y):
            memo[(n, m)] = True
            return True

    # If no move leads to a losing position for the opponent, this is a losing state.
    memo[(n, m)] = False
    return False

def find_p2_winning_starts(limit):
    """
    Finds the initial token counts (T) up to a given limit for which Player 2 has a winning strategy.
    """
    p2_wins = []
    for T in range(1, limit + 1):
        # Clear memoization cache for each new T to avoid excessively large cache.
        # This is a simple approach; more advanced caching could be used for larger limits.
        global memo
        memo = {}
        
        is_p1_win = False
        # Player 1 checks all possible first moves.
        for x1 in range(1, T):
            remaining_tokens = T - x1
            p2_max_move = 2 * x1
            
            # If P1 can make a move that leaves P2 in a losing position...
            if not is_winning_state(remaining_tokens, p2_max_move):
                # ...then T is a winning position for P1.
                is_p1_win = True
                break
        
        # If P1 could not find any such move, T is a winning position for P2.
        if not is_p1_win:
            p2_wins.append(T)
    return p2_wins

# We are looking for the T values where the second player wins.
# Let's find these values up to 100.
winning_Ts = find_p2_winning_starts(100)

print("The values of T for which the second player has a winning strategy (up to 100) are:")
for t in winning_Ts:
    print(t)