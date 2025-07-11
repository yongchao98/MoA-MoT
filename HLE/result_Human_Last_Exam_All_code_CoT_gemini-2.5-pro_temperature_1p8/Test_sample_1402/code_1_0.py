import math

def solve_theo_wins_probability():
    """
    Calculates the probability that Theo wins for the first time only after at least five games.

    The logic is as follows:
    1.  First, we determine the probability of Theo winning a single game, let's call it p_T.
        The game can be modeled as a 1D random walk on the difference d = (Heads - Tails).
        The walk starts at d=0. Absorbing states are d=+3 (Alex wins), d=-3 (Theo wins),
        and d=0 (Draw). By solving the standard gambler's ruin problem with these boundary
        conditions, we find that the probability of Alex winning, P(A), is 1/6.

    2.  By symmetry, the probability of Theo winning, p_T, is the same.
        p_T = P(A) = 1/6.

    3.  The probability that Theo does NOT win a single game is (1 - p_T).

    4.  The problem asks for the probability that Theo's first win occurs after at least five games.
        This is the same as the probability that Theo does not win in game 1, AND not in game 2,
        AND not in game 3, AND not in game 4.
        Since the games are independent events, we can multiply their probabilities.

    5.  The final probability is (1 - p_T)^4.
    """
    # Probability of Theo winning a single game
    p_T_win_num = 1
    p_T_win_den = 6

    # Probability of Theo not winning a single game
    p_T_not_win_num = p_T_win_den - p_T_win_num
    p_T_not_win_den = p_T_win_den

    # The event is "Theo does not win in the first 4 games"
    num_games = 4

    # Calculate the final probability as a fraction
    final_prob_num = p_T_not_win_num ** num_games
    final_prob_den = p_T_not_win_den ** num_games

    # Calculate the decimal value
    final_prob_decimal = final_prob_num / final_prob_den

    print(f"The probability of Theo winning a single game is {p_T_win_num}/{p_T_win_den}.")
    print(f"The probability of Theo not winning a single game is {p_T_not_win_num}/{p_T_not_win_den}.")
    print("The probability of Theo not winning for 4 consecutive games is calculated as:")
    print(f"Equation: ({p_T_not_win_num} / {p_T_not_win_den}) ^ {num_games} = {final_prob_num} / {final_prob_den}")
    print(f"\nThe final probability is {final_prob_num}/{final_prob_den}, which is approximately {final_prob_decimal:.8f}.")

solve_theo_wins_probability()
<<<625/1296>>>