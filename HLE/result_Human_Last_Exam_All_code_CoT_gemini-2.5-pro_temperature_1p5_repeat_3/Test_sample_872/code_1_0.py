from fractions import Fraction

def calculate_tic_tac_toe_win_probability():
    """
    Calculates the maximum win probability in Tic-Tac-Toe against a random opponent.

    The strategy involves:
    1. Player (X) starts in the center.
    2. Computer (O) plays randomly.
    3. We analyze two cases for the computer's first move: a corner or an edge.
    4. The total probability is the weighted average of the win probabilities of these cases.
    """

    # Probabilities of the computer's first move type, given you took the center.
    prob_O_plays_corner = Fraction(4, 8)
    prob_O_plays_edge = Fraction(4, 8)

    # Win probability if the computer plays a corner.
    # Calculation: (5/6) * 1 (computer doesn't block) + (1/6) * (7/8) (computer blocks, we continue optimally)
    win_prob_if_corner = Fraction(5, 6) * 1 + Fraction(1, 6) * Fraction(7, 8)

    # Win probability if the computer plays an edge.
    # Calculation: (5/6) * 1 (computer doesn't block) + (1/6) * 1 (computer blocks, we create a guaranteed fork)
    win_prob_if_edge = Fraction(5, 6) * 1 + Fraction(1, 6) * 1

    # Total win probability using the law of total probability.
    total_win_prob = (prob_O_plays_corner * win_prob_if_corner) + (prob_O_plays_edge * win_prob_if_edge)
    
    print("The optimal first move is to take the center square.")
    print("The total winning probability is calculated based on the computer's random response.")
    print("-" * 50)
    print("Final Winning Probability Calculation:")
    print(f"P(Win) = P(O plays corner) * P(Win|O plays corner) + P(O plays edge) * P(Win|O plays edge)")
    print(f"P(Win) = ({prob_O_plays_corner.numerator}/{prob_O_plays_corner.denominator}) * ({win_prob_if_corner.numerator}/{win_prob_if_corner.denominator}) + ({prob_O_plays_edge.numerator}/{prob_O_plays_edge.denominator}) * ({win_prob_if_edge.numerator}/{win_prob_if_edge.denominator})")
    
    term1 = prob_O_plays_corner * win_prob_if_corner
    term2 = prob_O_plays_edge * win_prob_if_edge
    
    print(f"P(Win) = ({term1.numerator}/{term1.denominator}) + ({term2.numerator}/{term2.denominator})")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")
    print("-" * 50)


calculate_tic_tac_toe_win_probability()