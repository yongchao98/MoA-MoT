import math

def solve_probability():
    """
    Calculates the probability that Theo wins for the first time only after at least five games.
    """
    
    # Step 1: Define the probabilities for a single game.
    # Through a state-based analysis of the game (a random walk), it can be shown that:
    # P(Theo wins) = 1/6
    # P(Alex wins) = 1/6
    # P(Draw) = 2/3
    # We are interested in the probability of Theo winning.
    p_theo_wins = 1/6
    
    # Step 2: The probability that Theo does NOT win a single game.
    p_theo_not_wins = 1 - p_theo_wins
    
    # The components of the fraction 5/6
    numerator = 5
    denominator = 6
    
    # Step 3: The problem asks for the probability that Theo's first win is on game 5 or later.
    # This is equivalent to Theo NOT winning the first four games.
    num_games = 4
    
    # Step 4: Calculate the final probability.
    # This is P(not win)^4 = (5/6)^4
    final_prob_num = numerator**num_games
    final_prob_den = denominator**num_games
    final_prob_decimal = final_prob_num / final_prob_den

    print("The probability of Theo winning a single game, P(T), is 1/6.")
    print("Therefore, the probability of Theo not winning a single game is 1 - 1/6 = 5/6.")
    print("\nThe event 'Theo wins for the first time only after at least five games' means Theo does not win the first four games.")
    print("The probability is calculated by raising the probability of not winning to the power of 4, since the games are independent.")
    
    print("\nThe final equation is:")
    print(f"P(win time >= 5) = ({numerator}/{denominator})^{num_games}")
    
    print("\nCalculating the result:")
    print(f"P = {final_prob_num} / {final_prob_den}")
    print(f"P ≈ {final_prob_decimal}")

solve_probability()