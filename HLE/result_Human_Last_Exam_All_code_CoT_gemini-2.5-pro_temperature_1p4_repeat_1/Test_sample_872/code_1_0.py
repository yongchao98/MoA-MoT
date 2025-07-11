from fractions import Fraction

def solve_tic_tac_toe_probability():
    """
    This function calculates and explains the maximum winning probability for Player X
    in Tic Tac Toe against a random opponent.
    """

    print("To find the maximum chance of winning, we must compare the outcomes of different opening moves.")
    print("We will analyze starting in a 'Corner' vs. starting in the 'Center'.\n")

    # --- Calculation for starting in a Corner ---
    print("--- Strategy 1: You Start in a Corner ---")
    print("When you start in a corner, the computer (Player O) has 8 possible moves.")
    print("Based on your optimal response to each, we find the following win probabilities:")
    print("- If O plays Center, Opposite Corner, or an Opposite Edge (4 of 8 moves): Your win is guaranteed. P(win) = 1")
    print("- If O plays an Adjacent Corner or an Adjacent Edge (4 of 8 moves): Your win probability is 23/24.\n")

    p_o_guaranteed_win = Fraction(4, 8)
    p_win_guaranteed = Fraction(1, 1)

    p_o_high_prob_win = Fraction(4, 8)
    p_win_high_prob = Fraction(23, 24)

    # Calculate total probability for starting in a corner
    P_corner = p_o_guaranteed_win * p_win_guaranteed + p_o_high_prob_win * p_win_high_prob
    
    print("The total win probability is calculated by summing the outcomes weighted by their likelihood:")
    print(f"P(win|corner) = (P(O's move leads to my guaranteed win) * 1) + (P(O's move leads to my high-prob win) * 23/24)")
    print(f"P(win|corner) = ({p_o_guaranteed_win.numerator}/{p_o_guaranteed_win.denominator}) * ({p_win_guaranteed.numerator}/{p_win_guaranteed.denominator}) + ({p_o_high_prob_win.numerator}/{p_o_high_prob_win.denominator}) * ({p_win_high_prob.numerator}/{p_win_high_prob.denominator})")
    print(f"P(win|corner) = {p_o_guaranteed_win} + {p_o_high_prob_win * p_win_high_prob}")
    print(f"P(win|corner) = {Fraction(1,2)} + {Fraction(23,48)}")
    print(f"P(win|corner) = {Fraction(24,48)} + {Fraction(23,48)} = {P_corner.numerator}/{P_corner.denominator}\n")

    # --- Calculation for starting in the Center ---
    print("--- Strategy 2: You Start in the Center ---")
    print("When you start in the center, the computer's 8 moves fall into two categories:")
    print("- If O plays a Corner (4 of 8 moves): Your win probability is 2/3.")
    print("- If O plays an Edge (4 of 8 moves): Your win is guaranteed. P(win) = 1.\n")

    p_o_corner_from_center = Fraction(4, 8)
    p_win_o_corner_from_center = Fraction(2, 3)

    p_o_edge_from_center = Fraction(4, 8)
    p_win_o_edge_from_center = Fraction(1, 1)

    # Calculate total probability for starting in the center
    P_center = p_o_corner_from_center * p_win_o_corner_from_center + p_o_edge_from_center * p_win_o_edge_from_center
    
    print("The total win probability is calculated as:")
    print(f"P(win|center) = (P(O plays corner) * 2/3) + (P(O plays edge) * 1)")
    print(f"P(win|center) = ({p_o_corner_from_center.numerator}/{p_o_corner_from_center.denominator}) * ({p_win_o_corner_from_center.numerator}/{p_win_o_corner_from_center.denominator}) + ({p_o_edge_from_center.numerator}/{p_o_edge_from_center.denominator}) * ({p_win_o_edge_from_center.numerator}/{p_win_o_edge_from_center.denominator})")
    print(f"P(win|center) = {p_o_corner_from_center * p_win_o_corner_from_center} + {p_o_edge_from_center * p_win_o_edge_from_center}")
    print(f"P(win|center) = {Fraction(1,3)} + {Fraction(1,2)} = {P_center.numerator}/{P_center.denominator}\n")

    # --- Comparison and Conclusion ---
    print("--- Conclusion ---")
    print(f"Comparing the probabilities: {P_corner.numerator}/{P_corner.denominator} (corner) vs. {P_center.numerator}/{P_center.denominator} (center).")
    
    # To compare, find a common denominator (48)
    P_center_common = P_center.limit_denominator(48)
    print(f"In other words: {P_corner.numerator}/{P_corner.denominator} vs. {P_center_common.numerator}/{P_center_common.denominator}.")
    
    max_prob = max(P_corner, P_center)
    print(f"\nThe maximum chance of winning you can give yourself is {max_prob.numerator}/{max_prob.denominator}.")

solve_tic_tac_toe_probability()
<<<47/48>>>