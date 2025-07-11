from fractions import Fraction

def solve_tic_tac_toe_probability():
    """
    Calculates the maximum win probability in Tic Tac Toe against a random opponent.

    The optimal first move is to take a corner. This function calculates the win
    probability by analyzing the computer's possible random responses and our
    optimal counter-moves.
    """
    print("Analyzing the optimal strategy: starting with X in a corner.")
    print("We calculate the total win probability by summing the weighted probabilities of winning after the computer's first move.\n")

    # Case A: Computer plays in the center.
    # Probability of this case: 1/8.
    # Our optimal play leads to a 5/6 chance of winning.
    prob_o1_center = Fraction(1, 8)
    win_prob_if_o1_center = Fraction(5, 6)
    win_contribution_A = prob_o1_center * win_prob_if_o1_center

    # Case B: Computer plays in the opposite corner.
    # Probability of this case: 1/8.
    # Our optimal play leads to a 2/3 chance of winning.
    prob_o1_opposite_corner = Fraction(1, 8)
    win_prob_if_o1_opposite_corner = Fraction(2, 3)
    win_contribution_B = prob_o1_opposite_corner * win_prob_if_o1_opposite_corner

    # Case C: Computer plays in an adjacent corner.
    # Probability of this case: 2/8.
    # Our optimal play guarantees a win.
    prob_o1_adjacent_corner = Fraction(2, 8)
    win_prob_if_o1_adjacent_corner = Fraction(1, 1)
    win_contribution_C = prob_o1_adjacent_corner * win_prob_if_o1_adjacent_corner

    # Case D: Computer plays in an adjacent edge.
    # Probability of this case: 2/8.
    # Our optimal play guarantees a win.
    prob_o1_adjacent_edge = Fraction(2, 8)
    win_prob_if_o1_adjacent_edge = Fraction(1, 1)
    win_contribution_D = prob_o1_adjacent_edge * win_prob_if_o1_adjacent_edge

    # Case E: Computer plays in a non-adjacent edge.
    # Probability of this case: 2/8.
    # Our optimal play guarantees a win.
    prob_o1_non_adjacent_edge = Fraction(2, 8)
    win_prob_if_o1_non_adjacent_edge = Fraction(1, 1)
    win_contribution_E = prob_o1_non_adjacent_edge * win_prob_if_o1_non_adjacent_edge
    
    print("The total win probability is the sum of probabilities of these cases:")
    print(f"P(Win) = P(O1=center) * P(Win|O1=center) +")
    print(f"         P(O1=opp corner) * P(Win|O1=opp corner) +")
    print(f"         P(O1=adj corner) * P(Win|O1=adj corner) +")
    print(f"         P(O1=adj edge) * P(Win|O1=adj edge) +")
    print(f"         P(O1=non-adj edge) * P(Win|O1=non-adj edge)\n")

    contributions = [win_contribution_A, win_contribution_B, win_contribution_C, win_contribution_D, win_contribution_E]
    total_win_prob = sum(contributions)

    # To show the common denominator for the equation
    common_denominator = 48
    parts_str = []
    for c in contributions:
        # Convert fraction to the common denominator for printing
        num = c.numerator * (common_denominator // c.denominator)
        parts_str.append(f"{num}/{common_denominator}")
    
    equation_str = " + ".join(parts_str)
    
    print("Substituting the fractional values into the equation:")
    print(f"P(Win) = ({prob_o1_center})*({win_prob_if_o1_center}) + ({prob_o1_opposite_corner})*({win_prob_if_o1_opposite_corner}) + ({prob_o1_adjacent_corner})*({win_prob_if_o1_adjacent_corner}) + ({prob_o1_adjacent_edge})*({win_prob_if_o1_adjacent_edge}) + ({prob_o1_non_adjacent_edge})*({win_prob_if_o1_non_adjacent_edge})")
    print(f"P(Win) = {win_contribution_A} + {win_contribution_B} + {win_contribution_C} + {win_contribution_D} + {win_contribution_E}")
    print(f"P(Win) = {equation_str}")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}\n")
    print(f"The maximum chance of winning is {total_win_prob}.")

solve_tic_tac_toe_probability()
<<<15/16>>>