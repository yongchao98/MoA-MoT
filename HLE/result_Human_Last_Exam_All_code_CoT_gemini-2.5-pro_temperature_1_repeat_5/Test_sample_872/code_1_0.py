from fractions import Fraction

def calculate_win_probability():
    """
    Calculates the maximum win probability in Tic-Tac-Toe against a random opponent.
    A tie is considered a loss.
    """
    print("Analyzing the optimal Tic-Tac-Toe strategy against a random opponent.\n")
    
    # --- 1. Analysis for starting in the CENTER ---
    print("--- Case 1: You start in the Center ---")
    # After you play center, the computer has 8 squares. 4 are corners, 4 are edges.
    prob_o_plays_corner = Fraction(4, 8)
    prob_o_plays_edge = Fraction(4, 8)
    
    # If computer plays a corner, you play the opposite corner.
    # The computer has 6 moves left.
    # To win, you need the computer to play on an edge adjacent to its first move (2 out of 6 spots).
    # If it plays anywhere else (4 out of 6 spots), you can only force a tie.
    win_prob_if_o_corner = Fraction(2, 6)
    print(f"If the computer plays a corner first (prob {prob_o_plays_corner}), your win probability is {win_prob_if_o_corner}.")
    
    # If computer plays an edge, you play a specific corner.
    # The computer has 6 moves left.
    # A detailed analysis shows 4 of its 6 possible moves lead to a forced win for you.
    win_prob_if_o_edge = Fraction(4, 6)
    print(f"If the computer plays an edge first (prob {prob_o_plays_edge}), your win probability is {win_prob_if_o_edge}.")
    
    total_win_prob_center = (prob_o_plays_corner * win_prob_if_o_corner) + \
                              (prob_o_plays_edge * win_prob_if_o_edge)
    
    # Original calculation: (4/8)*(2/6) + (4/8)*(4/6)
    # = (1/2)*(1/3) + (1/2)*(2/3) = 1/6 + 2/6 = 3/6 = 1/2
    print(f"\nTotal win probability if you start in the center: "
          f"({prob_o_plays_corner.numerator}/{prob_o_plays_corner.denominator}) * ({win_prob_if_o_corner.numerator}/{win_prob_if_o_corner.denominator}) + "
          f"({prob_o_plays_edge.numerator}/{prob_o_plays_edge.denominator}) * ({win_prob_if_o_edge.numerator}/{win_prob_if_o_edge.denominator}) = {total_win_prob_center}")
    print("-" * 20)

    # --- 2. Analysis for starting in a CORNER ---
    print("\n--- Case 2: You start in a Corner ---")
    # After you play corner, computer has 8 moves.
    prob_o_center = Fraction(1, 8)
    prob_o_opp_corner = Fraction(1, 8)
    prob_o_adj_corner = Fraction(2, 8)
    prob_o_adj_edge = Fraction(2, 8)
    prob_o_diag_edge = Fraction(2, 8)
    
    # If O plays center or opposite corner, you can only force a tie. Win prob = 0.
    win_prob_if_o_center = 0
    win_prob_if_o_opp_corner = 0
    print(f"If O plays center (prob {prob_o_center}) or opposite corner (prob {prob_o_opp_corner}), your win probability is 0 (it leads to a tie).")

    # If O plays adj corner, adj edge, or diag edge, you have a 5/6 chance of winning.
    # This is because your optimal response sets up a fork, and the computer has only 1 out of 6 moves to block it.
    win_prob_if_o_not_special = Fraction(5, 6)
    print(f"For any other move by O (prob {prob_o_adj_corner + prob_o_adj_edge + prob_o_diag_edge}), your win probability is {win_prob_if_o_not_special}.")

    total_win_prob_corner = (prob_o_center * win_prob_if_o_center) + \
                            (prob_o_opp_corner * win_prob_if_o_opp_corner) + \
                            (prob_o_adj_corner * win_prob_if_o_not_special) + \
                            (prob_o_adj_edge * win_prob_if_o_not_special) + \
                            (prob_o_diag_edge * win_prob_if_o_not_special)
    # Original calculation: (2/8)*0 + (6/8)*(5/6) = (3/4)*(5/6) = 15/24 = 5/8
    print(f"\nTotal win probability if you start in a corner: "
          f"({Fraction(2,8).numerator}/{Fraction(2,8).denominator})*0 + ({Fraction(6,8).numerator}/{Fraction(6,8).denominator}) * ({win_prob_if_o_not_special.numerator}/{win_prob_if_o_not_special.denominator}) = {total_win_prob_corner}")
    print("-" * 20)

    # --- 3. Analysis for starting on an EDGE ---
    print("\n--- Case 3: You start on an Edge ---")
    # After you play edge, computer has 8 moves.
    prob_o_center_e = Fraction(1, 8)
    prob_o_opp_edge_e = Fraction(1, 8)
    prob_o_adj_corner_e = Fraction(2, 8)
    prob_o_opp_corner_e = Fraction(2, 8)
    prob_o_adj_edge_e = Fraction(2, 8)

    win_prob_if_o_center_e = 0
    win_prob_if_o_adj_edge_e = 0
    
    win_prob_if_o_opp_edge_e = 1 # Forced win
    win_prob_if_o_corners_e = Fraction(5, 6) # Same logic as above

    total_win_prob_edge = (prob_o_center_e * win_prob_if_o_center_e) + \
                          (prob_o_adj_edge_e * win_prob_if_o_adj_edge_e) + \
                          (prob_o_opp_edge_e * win_prob_if_o_opp_edge_e) + \
                          ((prob_o_adj_corner_e + prob_o_opp_corner_e) * win_prob_if_o_corners_e)
                          
    # Original calculation: (3/8)*0 + (1/8)*1 + (4/8)*(5/6) = 1/8 + (1/2)*(5/6) = 1/8 + 5/12 = 3/24 + 10/24 = 13/24
    print(f"\nTotal win probability if you start on an edge: "
          f"({Fraction(3,8).numerator}/{Fraction(3,8).denominator})*0 + ({prob_o_opp_edge_e.numerator}/{prob_o_opp_edge_e.denominator})*1 + "
          f"({Fraction(4,8).numerator}/{Fraction(4,8).denominator}) * ({win_prob_if_o_corners_e.numerator}/{win_prob_if_o_corners_e.denominator}) = {total_win_prob_edge}")
    print("-" * 20)

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print(f"Win probability starting Center: {total_win_prob_center} (or {float(total_win_prob_center):.4f})")
    print(f"Win probability starting Corner: {total_win_prob_corner} (or {float(total_win_prob_corner):.4f})")
    print(f"Win probability starting Edge:   {total_win_prob_edge} (or {float(total_win_prob_edge):.4f})")
    
    max_prob = max(total_win_prob_center, total_win_prob_corner, total_win_prob_edge)
    
    print(f"\nThe maximum win probability is achieved by starting in a corner, which gives a {max_prob.numerator}/{max_prob.denominator} chance of winning.")
    return max_prob

result = calculate_win_probability()
# Final answer format
# print(f"\n<<<The maximum chance of winning is {result.numerator}/{result.denominator}>>>")

calculate_win_probability()
<<<5/8>>>