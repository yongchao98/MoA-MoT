def solve_coin_game_probability():
    """
    Calculates the probability that Theo wins for the first time
    only after at least five games.

    The solution is structured in three steps:
    1. Calculate the probability of Theo winning a single game (p_T).
    2. Calculate the probability of Theo not winning a single game.
    3. Calculate the probability of Theo not winning for the first four games.
    """

    print("--- Step 1: Probability of Theo winning a single game (p_T) ---")

    # For Theo to win, the sequence of tosses must result in T being 3 greater than H.
    # This can be modeled as a random walk on the state d = H - T.
    # The game starts at d=0. A head moves to d+1, a tail to d-1.
    # The game ends if d reaches 3 (Alex wins), -3 (Theo wins), or 0 (Draw).

    # The first toss must be a Tail (T) for Theo to have a chance to win.
    # If the first toss is a Head (H), the state becomes d=1. To reach d=-3, the walk
    # must cross d=0, which would end the game in a draw.
    prob_first_toss_T = 1 / 2

    # If the first toss is a Tail, the state is d=-1. Now, we need the probability
    # of the walk reaching d=-3 before it reaches d=0. This is a classic gambler's
    # ruin problem. The probability of reaching boundary 'b' before 'a' from state 'i'
    # is |i-a| / |b-a|. Here, i=-1, a=0, b=-3.
    prob_win_from_neg_1 = abs(-1 - 0) / abs(-3 - 0)

    # The total probability of Theo winning is the product of these two events.
    p_T = prob_first_toss_T * prob_win_from_neg_1
    
    print("For Theo to win, the first toss must be Tails (probability 1/2).")
    print("From the resulting state of -1, the probability of reaching -3 (Theo win) before 0 (Draw) is 1/3.")
    print(f"Thus, the probability of Theo winning a single game is p_T = (1/2) * (1/3) = {p_T:.4f} (or 1/6).")

    print("\n--- Step 2: Probability of Theo NOT winning a single game ---")
    prob_not_T = 1 - p_T
    # The probability of Theo not winning is simply 1 minus the probability of him winning.
    print(f"The probability that Theo does not win a game is 1 - p_T = 1 - 1/6 = 5/6.")

    print("\n--- Step 3: Probability over the first four games ---")
    # The event "Theo wins for the first time only after at least five games" means that
    # Theo did not win Game 1, Game 2, Game 3, AND Game 4.
    num_games_without_win = 4
    
    # Since each game is an independent event, we can multiply the probabilities.
    final_prob = prob_not_T ** num_games_without_win
    
    numerator = 5 ** num_games_without_win
    denominator = 6 ** num_games_without_win

    print("The event in question is equivalent to Theo not winning any of the first four games.")
    print("The probability for this is calculated as (P(Theo does not win))^4.")
    print("\n--- Final Calculation ---")
    print(f"The final equation is: ({5}/{6}) ^ {4} = {numerator}/{denominator}")
    print(f"The final probability as a decimal is approximately {final_prob:.6f}")

solve_coin_game_probability()