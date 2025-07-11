import fractions

def solve_tic_tac_toe_probability():
    """
    Calculates the maximum chance of winning Tic Tac Toe against a random computer.

    The analysis assumes we (Player X) go first, and a tie is considered a loss.
    Our strategy is to always make the move that maximizes our win probability.
    """

    print("Step 1: Determine the optimal opening move.")
    print("The optimal first move (X1) is to take the center square. This move is part of the most potential winning lines (4) and gives us the most control.\n")

    print("Step 2: Analyze the computer's possible first moves (O1).")
    print("After we take the center, the computer has 8 random choices. By symmetry, these can be grouped into two cases:")
    # Probabilities for Computer's First Move
    prob_o_corner = fractions.Fraction(4, 8)
    prob_o_edge = fractions.Fraction(4, 8)
    print(f"  - Case A: Computer plays in a corner (4 of 8 spots). Probability = {prob_o_corner.numerator}/{prob_o_corner.denominator}")
    print(f"  - Case B: Computer plays on an edge (4 of 8 spots). Probability = {prob_o_edge.numerator}/{prob_o_edge.denominator}\n")

    # --- Step 3a: Analysis of Case A (Computer plays a corner) ---
    print("--- Analysis of Case A (X=center, O=corner) ---")
    print("Our best response (X2) is to play on an edge. This creates an immediate threat of winning on our next turn.")
    print("The computer now has 6 empty squares for its second move (O2). To prevent our win, it must play in one specific square.")
    prob_o_fails_to_block_A = fractions.Fraction(5, 6)
    prob_o_blocks_A = fractions.Fraction(1, 6)
    
    # Sub-analysis for when the computer correctly blocks in Case A.
    # Our optimal play (X3) from there creates a complex situation for O3.
    # For O3's 4 random choices, 3 lead to a guaranteed win for us, and 1 leads to a 1/2 chance of winning.
    prob_win_if_o_blocks_A = fractions.Fraction(3, 4) * 1 + fractions.Fraction(1, 4) * fractions.Fraction(1, 2)
    
    print(f"The probability the computer FAILS to block is {prob_o_fails_to_block_A.numerator}/{prob_o_fails_to_block_A.denominator}. If it fails, we win.")
    print(f"The probability the computer BLOCKS is {prob_o_blocks_A.numerator}/{prob_o_blocks_A.denominator}.")
    print(f"If the computer blocks, our subsequent optimal play gives us a {prob_win_if_o_blocks_A.numerator}/{prob_win_if_o_blocks_A.denominator} chance of winning.")

    win_prob_A = prob_o_fails_to_block_A * 1 + prob_o_blocks_A * prob_win_if_o_blocks_A
    print(f"Total Win Probability for Case A = ({prob_o_fails_to_block_A}) + ({prob_o_blocks_A} * {prob_win_if_o_blocks_A}) = {win_prob_A}\n")

    # --- Step 3b: Analysis of Case B (Computer plays an edge) ---
    print("--- Analysis of Case B (X=center, O=edge) ---")
    print("Our best response (X2) is to play in a corner. This also creates an immediate threat.")
    print("The computer again has 6 empty squares and must block in one specific spot to prevent our win.")
    prob_o_fails_to_block_B = fractions.Fraction(5, 6)
    prob_o_blocks_B = fractions.Fraction(1, 6)
    
    # If the computer blocks, our next move (X3) creates a "fork" (two ways to win). A fork guarantees a win.
    prob_win_if_o_blocks_B = fractions.Fraction(1, 1)

    print(f"The probability the computer FAILS to block is {prob_o_fails_to_block_B.numerator}/{prob_o_fails_to_block_B.denominator}. If it fails, we win.")
    print(f"The probability the computer BLOCKS is {prob_o_blocks_B.numerator}/{prob_o_blocks_B.denominator}.")
    print(f"If the computer blocks, our next move creates a guaranteed win, making the win probability {prob_win_if_o_blocks_B}.")

    win_prob_B = prob_o_fails_to_block_B * 1 + prob_o_blocks_B * prob_win_if_o_blocks_B
    print(f"Total Win Probability for Case B = ({prob_o_fails_to_block_B}) + ({prob_o_blocks_B} * {prob_win_if_o_blocks_B}) = {win_prob_B}\n")

    # --- Step 4: Final Calculation ---
    print("--- Final Calculation ---")
    total_win_prob = prob_o_corner * win_prob_A + prob_o_edge * win_prob_B
    term_A = prob_o_corner * win_prob_A
    term_B = prob_o_edge * win_prob_B
    
    print("The total probability is the sum of probabilities of each case, weighted by their likelihood.")
    print("Total Win Probability = P(Case A) * P(Win in Case A) + P(Case B) * P(Win in Case B)")
    
    print("\nHere is the final equation with all the numbers:")
    print(f"({prob_o_corner.numerator}/{prob_o_corner.denominator}) * ({win_prob_A.numerator}/{win_prob_A.denominator}) + ({prob_o_edge.numerator}/{prob_o_edge.denominator}) * ({win_prob_B.numerator}/{win_prob_B.denominator}) "
          f"= ({term_A.numerator}/{term_A.denominator}) + ({term_B.numerator}/{term_B.denominator}) "
          f"= {total_win_prob.numerator}/{total_win_prob.denominator}")

    print("\nSo, the maximum chance of winning is given by the final fraction.")


solve_tic_tac_toe_probability()
<<<95/96>>>