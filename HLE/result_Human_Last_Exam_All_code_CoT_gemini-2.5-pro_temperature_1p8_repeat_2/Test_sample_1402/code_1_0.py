import math

def solve_coin_game_probability():
    """
    Calculates the probability that Theo wins for the first time only after at least five games.
    """

    # Step 1: Define the probabilities for a single game.
    # From analyzing the game as a state machine (or random walk), we can derive these probabilities.
    # P(Theo Wins) is the probability 'p' that Theo wins a single game.
    # The numerator and denominator for p are 1 and 6 respectively.
    p_theo_wins_num = 1
    p_theo_wins_den = 6

    # The probability that Theo does NOT win a single game is 1 - p.
    p_not_theo_win_num = p_theo_wins_den - p_theo_wins_num
    p_not_theo_win_den = p_theo_wins_den

    # Step 2: The problem asks for the probability that Theo wins for the first time only after
    # at least five games. This means Theo does not win the first four games.
    num_games = 4

    # Step 3: Calculate the final probability.
    # The outcomes of the games are independent, so we raise the probability of 'not winning'
    # to the power of the number of games.
    final_prob_numerator = p_not_theo_win_num ** num_games
    final_prob_denominator = p_not_theo_win_den ** num_games
    final_probability = final_prob_numerator / final_prob_denominator

    # Print the explanation and the final result.
    print("Let P(T) be the probability that Theo wins a single game.")
    print(f"Through analysis of the game's states, we find P(T) = {p_theo_wins_num}/{p_theo_wins_den}.")
    print("")
    print("The probability that Theo does NOT win a single game is 1 - P(T).")
    print(f"P(Not T) = 1 - {p_theo_wins_num}/{p_theo_wins_den} = {p_not_theo_win_num}/{p_not_theo_win_den}")
    print("")
    print("The event 'Theo wins for the first time after at least five games' is equivalent to 'Theo does not win any of the first four games'.")
    print(f"We need to calculate (P(Not T))^{num_games}.")
    print("")
    print("The equation for the final probability is:")
    print(f"({p_not_theo_win_num}/{p_not_theo_win_den})^{num_games} = {p_not_theo_win_num}^{num_games} / {p_not_theo_win_den}^{num_games}")
    print(f"Which evaluates to: {final_prob_numerator} / {final_prob_denominator}")
    print("")
    print(f"The final probability is {final_prob_numerator}/{final_prob_denominator}, which is approximately {final_probability:.5f}")

solve_coin_game_probability()