import sys
from fractions import Fraction

def solve_probability():
    """
    This function calculates the probability that Theo wins for the first time
    only after at least five games.
    """

    # Step 1: Calculate the probability of Theo winning a single game (p_T).
    # We model the game as a random walk on the states i = n_H - n_T.
    # The states are {-3, -2, -1, 0, 1, 2, 3}.
    # Theo wins at -3, Alex wins at 3, and it's a draw at 0.
    # Let Q(i) be the probability Theo wins starting from state i.
    # Boundary conditions: Q(-3) = 1, Q(3) = 0, Q(0) = 0.
    # For other states, Q(i) = 0.5 * Q(i-1) + 0.5 * Q(i+1).

    # We need to solve for Q(-1) and Q(1).
    # From Q(1) = 0.5*Q(0) + 0.5*Q(2) = 0.5*Q(2)
    # and Q(2) = 0.5*Q(1) + 0.5*Q(3) = 0.5*Q(1), we get Q(1)=0.
    #
    # From Q(-1) = 0.5*Q(-2) + 0.5*Q(0) = 0.5*Q(-2)
    # and Q(-2) = 0.5*Q(-3) + 0.5*Q(-1) = 0.5*1 + 0.5*Q(-1),
    # we can substitute: Q(-1) = 0.5 * (0.5 + 0.5*Q(-1)) => 0.75*Q(-1) = 0.25 => Q(-1) = 1/3.
    
    q_minus_1 = Fraction(1, 3)
    q_plus_1 = Fraction(0, 1)

    # The game starts at state 0. The first toss is either Heads (to state 1)
    # or Tails (to state -1), each with probability 1/2.
    # p_T = P(first toss T) * Q(-1) + P(first toss H) * Q(1)
    p_T = Fraction(1, 2) * q_minus_1 + Fraction(1, 2) * q_plus_1
    
    # Step 2: Calculate the probability of Theo NOT winning a single game.
    prob_not_win = 1 - p_T

    # Step 3: Calculate the probability of Theo not winning in the first 4 games.
    # The question asks for the probability that Theo's first win occurs on game 5 or later.
    # This is equivalent to Theo not winning games 1, 2, 3, and 4.
    num_games = 4
    final_prob = prob_not_win ** num_games
    
    # Step 4: Print the results clearly.
    print(f"The probability of Theo winning a single game is P(T) = {p_T.numerator}/{p_T.denominator}.")
    
    # Output the numbers in the final equation
    print("\nThe problem asks for the probability that Theo does not win in the first 4 games.")
    print(f"The probability of Theo NOT winning a single game is 1 - P(T) = {prob_not_win.numerator}/{prob_not_win.denominator}.")
    
    print("\nThe final equation is (numerator / denominator) ** exponent:")
    print(f"  Numerator: {prob_not_win.numerator}")
    print(f"  Denominator: {prob_not_win.denominator}")
    print(f"  Exponent: {num_games}")

    # Final Result
    print(f"\nThe final probability is ({prob_not_win.numerator}/{prob_not_win.denominator})^{num_games} = {final_prob.numerator}/{final_prob.denominator}.")
    print(f"As a decimal, this is approximately {float(final_prob):.6f}.")

solve_probability()