from fractions import Fraction

def calculate_tic_tac_toe_win_probability():
    """
    Calculates and explains the maximum probability of winning Tic Tac Toe
    against a computer that plays randomly.
    """
    # Player X starts by taking the center square. This is the optimal first move.
    # The computer (Player O) then has 8 empty squares to choose from.
    # By symmetry, there are two cases for the computer's first move:
    # 1. A corner square (4 available corners out of 8 squares).
    # 2. An edge square (4 available edges out of 8 squares).

    prob_O_chooses_corner = Fraction(4, 8)
    prob_O_chooses_edge = Fraction(4, 8)

    # --- Case 1: Computer chooses a corner square ---
    # After X's optimal response to create a threat, there are 6 empty squares.
    # O must place a mark in the single square that blocks X's threat.
    prob_O_fails_first_block = Fraction(5, 6)
    prob_O_succeeds_first_block = Fraction(1, 6)

    # If O succeeds, X makes another threat. There are now 4 empty squares.
    prob_O_fails_second_block = Fraction(3, 4)

    # The win probability is P(O fails first) + P(O succeeds first) * P(O fails second).
    # If O blocks both, the game is a draw (not a win).
    win_prob_if_O_takes_corner = prob_O_fails_first_block + prob_O_succeeds_first_block * prob_O_fails_second_block

    # --- Case 2: Computer chooses an edge square ---
    # X's optimal move creates a threat that, if blocked, sets up a guaranteed win (a "fork").
    # After X's move, O has 6 squares left and must block the immediate threat.
    prob_O_fails_to_block_setup = Fraction(5, 6)
    prob_O_succeeds_in_blocking_setup = Fraction(1, 6)

    # If O blocks, X's next move creates a fork, guaranteeing a win.
    win_prob_after_block = Fraction(1, 1)
    win_prob_if_O_takes_edge = prob_O_fails_to_block_setup + prob_O_succeeds_in_blocking_setup * win_prob_after_block

    # --- Total Probability ---
    total_win_prob = (prob_O_chooses_corner * win_prob_if_O_takes_corner) + (prob_O_chooses_edge * win_prob_if_O_takes_edge)

    # --- Output the calculation and result ---
    print("The maximum chance of winning is calculated by starting in the center.")
    print("The calculation depends on the computer's random response (playing a corner or an edge).\n")
    print("P(Win) = P(O plays corner) * P(Win | O plays corner) + P(O plays edge) * P(Win | O plays edge)\n")
    
    # Print the equation with all the numbers
    print("The final equation is:")
    print(f"P(Win) = ({prob_O_chooses_corner.numerator}/{prob_O_chooses_corner.denominator}) * ({prob_O_fails_first_block.numerator}/{prob_O_fails_first_block.denominator} + ({prob_O_succeeds_first_block.numerator}/{prob_O_succeeds_first_block.denominator} * {prob_O_fails_second_block.numerator}/{prob_O_fails_second_block.denominator})) + ({prob_O_chooses_edge.numerator}/{prob_O_chooses_edge.denominator}) * ({prob_O_fails_to_block_setup.numerator}/{prob_O_fails_to_block_setup.denominator} + ({prob_O_succeeds_in_blocking_setup.numerator}/{prob_O_succeeds_in_blocking_setup.denominator} * {win_prob_after_block.numerator}/{win_prob_after_block.denominator}))")
    print(f"P(Win) = ({prob_O_chooses_corner.numerator}/{prob_O_chooses_corner.denominator}) * ({win_prob_if_O_takes_corner.numerator}/{win_prob_if_O_takes_corner.denominator}) + ({prob_O_chooses_edge.numerator}/{prob_O_chooses_edge.denominator}) * ({win_prob_if_O_takes_edge.numerator}/{win_prob_if_O_takes_edge.denominator})")
    
    term1 = prob_O_chooses_corner * win_prob_if_O_takes_corner
    term2 = prob_O_chooses_edge * win_prob_if_O_takes_edge
    
    print(f"P(Win) = {term1.numerator}/{term1.denominator} + {term2.numerator}/{term2.denominator}")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}\n")
    print(f"The maximum chance of winning is {total_win_prob.numerator}/{total_win_prob.denominator}.")

calculate_tic_tac_toe_win_probability()