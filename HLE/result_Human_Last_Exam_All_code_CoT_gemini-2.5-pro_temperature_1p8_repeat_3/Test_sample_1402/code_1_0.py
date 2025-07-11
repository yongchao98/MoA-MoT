import math

def solve_probability():
    """
    Calculates the probability that Theo wins for the first time only after at least five games.
    """
    # Step 1: Determine the probabilities for a single game outcome.
    # We model the game as a random walk on the state D = (Number of Heads - Number of Tails).
    # The game stops if D reaches 3 (Alex wins), -3 (Theo wins), or 0 (Draw).
    # Let p_i be the probability Alex wins starting from state D=i.
    # The game starts at D=0 and moves to D=1 or D=-1 after the first toss.
    # We can set up a system of linear equations:
    # 1) p_1 = 0.5 * p_2 + 0.5 * p_0_win
    #    When the state returns to 0, it's a draw, so Alex's win probability from that path is 0.
    #    So, p_1 = 0.5 * p_2
    # 2) p_2 = 0.5 * p_3 + 0.5 * p_1
    #    When the state reaches 3, Alex wins, so the probability of winning is 1.
    #    So, p_2 = 0.5 * 1 + 0.5 * p_1
    #
    # Solving for p_1:
    # Substitute (1) into (2): p_2 = 0.5 + 0.5 * (0.5 * p_2) => p_2 = 0.5 + 0.25 * p_2
    # 0.75 * p_2 = 0.5 => p_2 = 0.5 / 0.75 = 2/3
    # p_1 = 0.5 * p_2 = 0.5 * (2/3) = 1/3
    p_a_from_state_1 = 1/3

    # The overall probability of Alex winning, P_A, is P(first toss is H) * p_1
    prob_alex_wins = 0.5 * p_a_from_state_1

    # By symmetry, the probability of Theo winning is the same.
    prob_theo_wins = prob_alex_wins

    print("Step 1: Calculating the probabilities for a single game.")
    print(f"The probability of Alex winning a single game is 1/6, or approximately {prob_alex_wins:.4f}.")
    print(f"By symmetry, the probability of Theo winning a single game is also 1/6, or approximately {prob_theo_wins:.4f}.")
    
    # Step 2: Formulate the probability for the series of games.
    # The question is the probability that Theo wins for the first time on game k, where k >= 5.
    # This is equivalent to Theo NOT winning in the first 4 games.
    prob_theo_not_wins = 1 - prob_theo_wins
    
    # The number of games Theo must not win consecutively at the start.
    num_games = 4
    
    # The probability is P(Theo not wins)^4, as the games are independent.
    final_prob = prob_theo_not_wins ** num_games
    
    print("\nStep 2: Calculating the probability over multiple games.")
    print("The event 'Theo wins for the first time only after at least five games' means Theo must not win in the first four games.")
    print(f"The probability of Theo NOT winning a single game is 1 - 1/6 = 5/6.")
    print(f"Since the games are independent, the total probability is calculated by raising this to the power of {num_games}.")

    # Step 3: Display the final equation and result.
    base_numerator = 5
    base_denominator = 6
    exponent = 4
    final_numerator = base_numerator ** exponent
    final_denominator = base_denominator ** exponent
    
    print("\n--- Final Answer ---")
    print(f"The final probability is calculated with the equation: P = ({base_numerator}/{base_denominator})^{exponent}")
    print(f"The numbers in this equation are the base numerator {base_numerator}, the base denominator {base_denominator}, and the exponent {exponent}.")
    print(f"This evaluates to the fraction: {final_numerator}/{final_denominator}")
    print(f"The decimal value is approximately: {final_prob:.7f}")

solve_probability()