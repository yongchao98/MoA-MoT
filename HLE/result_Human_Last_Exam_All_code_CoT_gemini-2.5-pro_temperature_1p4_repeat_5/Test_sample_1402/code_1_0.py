import math

def solve_probability():
    """
    Calculates the probability that Theo wins for the first time
    only after at least five games.
    """

    # Step 1: Define the probabilities for a single game.
    # p_T: Probability Theo wins a single game. Based on the random walk analysis, p_T = 1/6.
    p_T_num = 1
    p_T_den = 6

    # The probability of Theo NOT winning a single game is 1 - p_T.
    p_not_T_num = p_T_den - p_T_num
    p_not_T_den = p_T_den

    # Step 2: Calculate the probability over a sequence of games.
    # The event "Theo wins for the first time only after at least five games"
    # is equivalent to "Theo does not win in any of the first four games".
    num_games = 4

    # The probability of this sequence is (1 - p_T)^4 = (5/6)^4.
    result_num = p_not_T_num ** num_games
    result_den = p_not_T_den ** num_games

    # Print the explanation and the final result.
    print("The problem asks for the probability that Theo's first win occurs at game k, where k >= 5.")
    print("This is equivalent to the probability that Theo does not win in games 1, 2, 3, and 4.")
    print(f"The probability of Theo winning a single game is {p_T_num}/{p_T_den}.")
    print(f"The probability of Theo not winning a single game is {p_not_T_num}/{p_not_T_den}.")
    print(f"Since each game is an independent event, the total probability is calculated by multiplying the probabilities for the first {num_games} games.")
    print("\nFinal Equation:")
    print(f"({p_not_T_num}/{p_not_T_den})^{num_games} = {result_num}/{result_den}")
    print(f"\nThe decimal value of the probability is approximately: {result_num / result_den:.6f}")

solve_probability()