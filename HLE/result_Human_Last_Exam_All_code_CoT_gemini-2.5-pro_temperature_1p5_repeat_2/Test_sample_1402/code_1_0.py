import fractions

def solve_coin_game_probability():
    """
    Calculates the probability that Theo wins for the first time only after at least five games.
    """

    # --- Step 1: Analyze a single game to find the probability of Theo winning (p_theo) ---

    # The game can be modeled as a random walk on the integer line, where the position
    # is the difference `d = (number of heads) - (number of tails)`.
    # The game starts at d=0. Alex wins if d=3, Theo wins if d=-3, and it's a draw if d returns to 0.

    # For Theo to win, the sequence of states must reach -3 without returning to 0 first.
    # This requires the first toss to be a Tail (T), moving the state from 0 to -1.
    prob_first_toss_T = fractions.Fraction(1, 2)

    # From state d=-1, we need to find the probability of reaching d=-3 before reaching d=0.
    # This is a classic "Gambler's Ruin" problem. For a fair game (p=0.5), a player starting
    # with fortune `i` out of a total `N` between two absorbing barriers has a probability `i/N`
    # of reaching the barrier at N.
    # Here, the 'fortune' is the distance from the 'draw' barrier (0). We start at d=-1, so our
    # initial fortune is i=1. The 'win' barrier for Theo is at d=-3, so the total distance is N=3.
    # The probability of reaching -3 from -1 is i/N.
    prob_win_from_neg_1 = fractions.Fraction(1, 3)

    # The total probability of Theo winning a single game is the product of these probabilities.
    p_theo = prob_first_toss_T * prob_win_from_neg_1

    # --- Step 2: Calculate the probability over multiple games ---

    # The problem asks for the probability that Theo wins for the first time *after at least five games*.
    # This is equivalent to Theo NOT winning in the first four games.
    # The outcome of game 5 and beyond doesn't change the fact that the first four were not wins for Theo.

    # Probability of Theo NOT winning a single game.
    p_not_theo = 1 - p_theo

    # The number of initial games Theo must not win.
    num_games = 4

    # Since each game is an independent event, we raise the probability of not winning to the power of 4.
    final_prob = p_not_theo ** num_games

    # --- Step 3: Print the logic and the final result ---

    print("Step 1: Probability of Theo winning a single game (p_theo)")
    print(f"p_theo = P(first toss is Tails) * P(reach -3 from -1 before 0)")
    print(f"p_theo = (1/2) * (1/3) = {p_theo.numerator}/{p_theo.denominator}\n")

    print("Step 2: Probability of Theo not winning in the first 4 games")
    print("The event 'Theo wins for the first time after at least five games' means Theo does not win in games 1, 2, 3, and 4.")
    print(f"The probability of Theo NOT winning one game is 1 - {p_theo.numerator}/{p_theo.denominator} = {p_not_theo.numerator}/{p_not_theo.denominator}\n")

    print(f"The total probability is the probability of this happening {num_games} times in a row.")
    print("Final Equation:")
    print(f"({p_not_theo.numerator}/{p_not_theo.denominator})^{num_games} = "
          f"{p_not_theo.numerator}^{num_games} / {p_not_theo.denominator}^{num_games} = "
          f"{final_prob.numerator} / {final_prob.denominator}")

    print(f"\nThe decimal equivalent is approximately {float(final_prob):.6f}")

solve_coin_game_probability()