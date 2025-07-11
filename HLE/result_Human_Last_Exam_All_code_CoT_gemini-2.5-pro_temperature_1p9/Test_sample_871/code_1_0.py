import sys

# It's possible for N to be large enough to hit the default recursion limit.
# We can increase it to be safe, though it may not be necessary for this problem.
sys.setrecursionlimit(2000)

memo = {}
target_N = 0

def can_win(current_sum, last_move):
    """
    Determines if the current player can force a win from the given state.
    A state is defined by (current_sum, last_move).
    Memoization is used to store results for previously seen states to avoid re-computation.
    """
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # The current player can choose any number from 1 to 9, except the one
    # the opponent just played.
    possible_moves = list(range(1, 10))
    if last_move is not None and last_move in possible_moves:
        possible_moves.remove(last_move)
    
    # A position is a "winning" position if there is at least one move
    # that leads to a "losing" position for the opponent.
    for move in possible_moves:
        new_sum = current_sum + move
        
        # If a move reaches the target sum exactly, it's a winning move.
        if new_sum == target_N:
            memo[(current_sum, last_move)] = True
            return True
        
        if new_sum < target_N:
            # Check if this move forces the opponent into a losing position.
            # The opponent will be at state (new_sum, move). If they cannot win
            # from there, this is a good move for us.
            if not can_win(new_sum, move):
                memo[(current_sum, last_move)] = True
                return True

    # If no move can lead to a win (either they overshoot the target or
    # they all lead to winning positions for the opponent), this is a "losing" position.
    memo[(current_sum, last_move)] = False
    return False

def check_b_win_for_n(N):
    """
    Checks if Player B has a guaranteed winning strategy for a target sum N.
    """
    global memo, target_N
    memo.clear()
    target_N = N
    
    # Player B wins if, for every possible first move by Player A,
    # Player B can force a win from the resulting state.
    for first_move in range(1, 10):
        # After Player A's move 'first_move', the sum is 'first_move',
        # and it becomes Player B's turn. The last move was 'first_move'.
        # If Player B CANNOT win from this state, then Player A has found a
        # winning strategy, so this N is not a win for B.
        if not can_win(first_move, first_move):
            return False
            
    # If Player B can win regardless of Player A's first move, it's a B-win N.
    return True

def find_solution():
    """
    Finds the first three values of N for which Player B wins and prints their product.
    """
    b_winning_ns = []
    n = 1
    while len(b_winning_ns) < 3:
        if check_b_win_for_n(n):
            b_winning_ns.append(n)
        n += 1
        
    product = 1
    for val in b_winning_ns:
        product *= val
        
    n1, n2, n3 = b_winning_ns
    print(f"The three smallest values of N for which Player B can win are: {n1}, {n2}, {n3}")
    print(f"The final equation is: {n1} * {n2} * {n3} = {product}")
    # The final answer in the required format.
    print(f"\n<<<{product}>>>")

if __name__ == '__main__':
    find_solution()
