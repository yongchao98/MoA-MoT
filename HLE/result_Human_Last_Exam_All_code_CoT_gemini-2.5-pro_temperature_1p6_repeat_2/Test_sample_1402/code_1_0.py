import math

def solve_theo_win_probability():
    """
    Calculates the probability that Theo wins for the first time only after at least five games.

    The steps are:
    1. Determine the probability of Theo winning a single game (p_T).
       Based on a random walk model, this probability is 1/6.
    2. Determine the probability of Theo NOT winning a single game (1 - p_T).
       This is 1 - 1/6 = 5/6.
    3. The event "Theo wins for the first time after at least five games" is equivalent
       to "Theo does not win in game 1, 2, 3, and 4".
    4. Calculate the probability of this sequence of independent events: (5/6)^4.
    """

    # Probability of Theo winning a single game, as derived from the random walk model.
    p_T_numerator = 1
    p_T_denominator = 6

    # Probability of Theo NOT winning a single game.
    p_not_T_numerator = p_T_denominator - p_T_numerator
    p_not_T_denominator = p_T_denominator

    # Number of games Theo must not win consecutively.
    num_games = 4

    # Calculate the numerator and denominator of the final probability, which is (5/6)^4.
    final_numerator = p_not_T_numerator ** num_games
    final_denominator = p_not_T_denominator ** num_games

    # Calculate the decimal value of the result.
    final_probability = final_numerator / final_denominator

    print("This script calculates the probability that Theo's first win occurs after at least five games.")
    print("-" * 80)
    print(f"The probability of Theo winning a single game is {p_T_numerator}/{p_T_denominator}.")
    print(f"The probability of Theo NOT winning a single game is {p_not_T_numerator}/{p_not_T_denominator}.")
    print(f"The number of games Theo must not win for the condition to hold is {num_games}.")
    print("\nThe final probability is the probability of not winning raised to the power of the number of games.")
    print("\nFinal Equation:")
    print(f"({p_not_T_numerator}/{p_not_T_denominator})^{num_games} = ({p_not_T_numerator}^{num_games}) / ({p_not_T_denominator}^{num_games}) = {final_numerator} / {final_denominator}")
    print("\nFinal Answer:")
    print(f"The probability is {final_numerator}/{final_denominator}, which is approximately {final_probability:.7f}")

solve_theo_win_probability()

# The final result as a decimal value for submission.
final_answer = (5/6)**4
# <<<...>>> format is for the final answer, not part of the code execution.
# So the code will not print the final answer tag.