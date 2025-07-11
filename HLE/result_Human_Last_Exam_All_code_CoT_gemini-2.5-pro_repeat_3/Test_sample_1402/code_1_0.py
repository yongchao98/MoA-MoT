import fractions

def solve_probability():
    """
    This function calculates the probability that Theo wins for the first time
    only after at least five games.
    """
    
    # Step 1: Calculate the probability of Theo winning a single game (P_T).
    #
    # We model the game as a random walk on the state d = H - T. The game ends
    # when d reaches 3 (Alex wins), -3 (Theo wins), or 0 (Draw).
    #
    # Let l_i be the probability Theo wins (reaches -3) starting from state d=i.
    # The boundary conditions are l_{-3} = 1, l_0 = 0, and l_3 = 0.
    # We have the recurrence relations for this 1D random walk:
    # l_i = 0.5 * l_{i+1} + 0.5 * l_{i-1}
    #
    # Solving for i < 0 (states -1, -2):
    # l_{-1} = 0.5 * l_0 + 0.5 * l_{-2} = 0.5 * l_{-2}
    # l_{-2} = 0.5 * l_{-1} + 0.5 * l_{-3} = 0.5 * l_{-1} + 0.5 * 1
    # Substituting l_{-1} into the second equation:
    # l_{-2} = 0.5 * (0.5 * l_{-2}) + 0.5 = 0.25 * l_{-2} + 0.5
    # 0.75 * l_{-2} = 0.5  => l_{-2} = 0.5 / 0.75 = 2/3
    # l_{-1} = 0.5 * l_{-2} = 0.5 * (2/3) = 1/3
    #
    # For i > 0, any path to -3 must cross 0, which is an absorbing state (draw).
    # So, l_1 = 0 and l_2 = 0.
    #
    # The game starts at d=0, and the first toss leads to d=1 (Heads) or d=-1 (Tails),
    # each with probability 0.5.
    # The total probability of Theo winning is:
    # P_T = 0.5 * l_1 + 0.5 * l_{-1} = 0.5 * 0 + 0.5 * (1/3) = 1/6
    
    prob_theo_wins_one_game = fractions.Fraction(1, 6)
    
    # Step 2: Calculate the probability over a sequence of games.
    # The event is "Theo wins for the first time only after at least five games".
    # This means Theo does not win in game 1, AND not in game 2, AND not in game 3,
    # AND not in game 4. The outcomes of the games are independent.
    
    # The probability that Theo does NOT win a single game is 1 - P_T.
    prob_theo_not_wins_one_game = 1 - prob_theo_wins_one_game
    
    # The number of initial games Theo must not win is 4.
    num_games_without_win = 4
    
    # The final probability is (P(Theo does not win))^4.
    final_prob = prob_theo_not_wins_one_game ** num_games_without_win
    
    # Print the explanation and the final equation.
    print("The probability of Theo winning a single game is 1/6.")
    print("Therefore, the probability of Theo not winning a single game is 1 - 1/6 = 5/6.")
    print("\nThe event 'Theo wins for the first time only after at least five games' is equivalent to 'Theo does not win in the first four games'.")
    print("Since each game is an independent event, we calculate this probability as (5/6)^4.")
    
    # Output the numbers in the final equation.
    numerator = prob_theo_not_wins_one_game.numerator
    denominator = prob_theo_not_wins_one_game.denominator
    power = num_games_without_win
    result_numerator = final_prob.numerator
    result_denominator = final_prob.denominator
    
    print("\nFinal Calculation:")
    print(f"({numerator}/{denominator})^{power} = {numerator**power}/{denominator**power} = {result_numerator}/{result_denominator}")

solve_probability()
<<<625/1296>>>