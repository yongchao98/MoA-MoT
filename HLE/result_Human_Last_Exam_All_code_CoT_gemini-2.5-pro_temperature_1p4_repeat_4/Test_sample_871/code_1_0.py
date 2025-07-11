import sys

# Set a higher recursion limit for safety, although the required N values are small enough
# that it's likely not necessary.
sys.setrecursionlimit(2000)

# Memoization table to store results of subproblems
memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can force a win from the given state (current_sum, last_move).
    Uses memoization to avoid re-computing results for the same state.
    
    A player can win if they can make a move that leads to a state from which
    the opponent cannot win.
    """
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # A move is any integer from 1 to 9 not equal to the opponent's last move.
    for k in range(1, 10):
        if k == last_move:
            continue
        
        new_sum = current_sum + k

        if new_sum == N:
            # This move immediately wins the game.
            memo[(current_sum, last_move)] = True
            return True
        
        if new_sum < N:
            # If we make this move, the opponent will be in the state (new_sum, k).
            # If the opponent CANNOT win from that new state, it means this move
            # leads to a losing position for them, so it's a winning move for us.
            if not can_win(new_sum, k, N):
                memo[(current_sum, last_move)] = True
                return True

    # If we've looped through all possible moves and none of them lead to a
    # losing position for the opponent, then the current position is a losing one for us.
    # This also handles the case where any move k makes new_sum > N (a loss).
    memo[(current_sum, last_move)] = False
    return False

def solve():
    """
    Finds the three smallest values of N for which Player B can win,
    then computes and prints their product.
    """
    b_win_Ns = []
    n = 1
    while len(b_win_Ns) < 3:
        # Player B wins if Player A is in a losing position at the start.
        # Player A starts at (sum=0, last_move=0).
        # So, we check if can_win(0, 0, n) is False for Player A.
        global memo
        memo = {}
        if not can_win(0, 0, n):
            b_win_Ns.append(n)
        n += 1

    product = b_win_Ns[0] * b_win_Ns[1] * b_win_Ns[2]
    
    print(f"The three smallest values of N for which Player B has a winning strategy are: {b_win_Ns[0]}, {b_win_Ns[1]}, and {b_win_Ns[2]}.")
    print("The product is:")
    print(f"{b_win_Ns[0]} * {b_win_Ns[1]} * {b_win_Ns[2]} = {product}")

solve()
<<<7986>>>