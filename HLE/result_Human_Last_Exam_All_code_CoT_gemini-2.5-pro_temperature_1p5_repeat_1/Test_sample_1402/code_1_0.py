import math

def solve_probability():
    """
    This function calculates the probability that Theo wins for the first time only after at least five games.
    """

    # Step 1: Define the probability of Theo winning a single game.
    # As explained in the plan, the probability of Theo winning (P_T), Alex winning (P_A),
    # and a draw (P_D) are:
    # P_A = 1/6
    # P_T = 1/6
    # P_D = 2/3
    p_theo_wins_num = 1
    p_theo_wins_den = 6
    
    # Step 2: Calculate the probability of Theo NOT winning a single game.
    # P(Theo does not win) = 1 - P(Theo wins) = 1 - 1/6 = 5/6
    p_theo_not_wins_num = p_theo_wins_den - p_theo_wins_num
    p_theo_not_wins_den = p_theo_wins_den
    
    print(f"The probability of Theo winning a single game is {p_theo_wins_num}/{p_theo_wins_den}.")
    print(f"The probability of Theo not winning a single game is {p_theo_not_wins_num}/{p_theo_not_wins_den}.")
    print("-" * 30)

    # Step 3: The event "Theo wins for the first time after at least five games" means
    # Theo does not win in game 1, AND not in game 2, AND not in game 3, AND not in game 4.
    # The probability is P(Theo does not win)^4 = (5/6)^4.
    num_games = 4
    
    # Step 4: Calculate the numerator and denominator of the final probability.
    final_prob_num = p_theo_not_wins_num ** num_games
    final_prob_den = p_theo_not_wins_den ** num_games
    
    # Step 5: Print the final equation and the result.
    print(f"The probability of Theo not winning in the first {num_games} games is calculated as ({p_theo_not_wins_num}/{p_theo_not_wins_den})^{num_games}.")
    
    print("\nThe final equation is:")
    # Using f-string to output each number in the equation
    print(f"{p_theo_not_wins_num}^{num_games} / {p_theo_not_wins_den}^{num_games} = {final_prob_num} / {final_prob_den}")

    final_probability = final_prob_num / final_prob_den
    print(f"\nThe probability as a decimal is approximately: {final_probability:.5f}")


solve_probability()
<<<625/1296>>>