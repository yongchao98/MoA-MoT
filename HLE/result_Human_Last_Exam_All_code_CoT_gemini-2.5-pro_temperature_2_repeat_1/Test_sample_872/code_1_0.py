from fractions import Fraction

def solve_tic_tac_toe_probability():
    """
    Calculates the maximum win probability for a player against a random-move computer
    in Tic Tac Toe, where the player goes first and a tie is considered a loss.
    """

    # --- Step 1: Analyze Computer's First Move ---
    # After you (X) take the center, there are 8 squares left for the computer (O).
    # 4 are corners and 4 are edges.
    p_O1_is_corner = Fraction(4, 8)
    p_O1_is_edge = Fraction(4, 8)

    # --- Step 2: Calculate Win Probability if Computer picks an EDGE first ---
    # Your best response is to take any corner. This creates a winning threat.
    # The computer has 6 remaining squares for its next move (O2).
    # It must block your threat to prevent an immediate loss.
    # The chance it FAILS to block (and you win immediately) is 5/6.
    p_O2_fails_to_block = Fraction(5, 6)
    
    # If the computer DOES block (1/6 chance), you can make a move that creates a "fork"
    # (two simultaneous winning threats). The random computer can only block one,
    # so you are guaranteed to win on the next turn.
    p_win_after_block = Fraction(1, 1)

    # Your win probability in this branch is:
    # P(Win | O1 is Edge) = P(O fails to block) * 1 + P(O blocks) * P(Win after block)
    #                   = 5/6 * 1 + 1/6 * 1 = 1
    p_win_if_O1_edge = Fraction(1, 1)

    # --- Step 3: Calculate Win Probability if Computer picks a CORNER first ---
    # Your best response is the opposite corner. This sets you up for a fork.
    # The computer has 6 squares for its move (O2). To stop your fork, it must
    # play on one of two specific squares.
    # The chance it AVOIDS those squares (letting you win) is 4/6.
    p_O2_avoids_fork_setup = Fraction(4, 6)

    # If it DOES land on a fork-blocking square (2/6 chance), you are forced to block its threat.
    # This move also creates a new threat for you. The computer has 4 squares for its
    # move (O3) and must block your new threat.
    # The chance it FAILS to block is 3/4.
    p_O3_fails_to_block_2nd_threat = Fraction(3, 4)

    # If it DOES block your 2nd threat (1/4 chance), you make a new threat.
    # The computer has 2 squares for its move (O4) and must block.
    # The chance it FAILS to block is 1/2. If it succeeds, the game is a tie (loss for you).
    p_O4_fails_to_block_3rd_threat = Fraction(1, 2)
    p_win_after_2nd_block = p_O4_fails_to_block_3rd_threat * 1 + (1 - p_O4_fails_to_block_3rd_threat) * 0

    # Combining these gives the win probability if your initial fork setup was blocked.
    p_win_after_fork_setup_blocked = p_O3_fails_to_block_2nd_threat * 1 + (1 - p_O3_fails_to_block_2nd_threat) * p_win_after_2nd_block
    
    # Your win probability in the corner branch is:
    p_win_if_O1_corner = p_O2_avoids_fork_setup * 1 + (1 - p_O2_avoids_fork_setup) * p_win_after_fork_setup_blocked

    # --- Step 4: Final Calculation ---
    # P(Win) = P(O1 is Edge) * P(Win | O1 is Edge) + P(O1 is Corner) * P(Win | O1 is Corner)
    total_win_prob = p_O1_is_edge * p_win_if_O1_edge + p_O1_is_corner * p_win_if_O1_corner

    print("The maximum chance of winning is calculated as follows:")
    print("P(Win) = P(Computer's 1st move is Edge) * P(Win | O1 is Edge) + P(Computer's 1st move is Corner) * P(Win | O1 is Corner)")
    print(f"P(Win) = {p_O1_is_edge.numerator}/{p_O1_is_edge.denominator} * {p_win_if_O1_edge.numerator}/{p_win_if_O1_edge.denominator} + {p_O1_is_corner.numerator}/{p_O1_is_corner.denominator} * {p_win_if_O1_corner.numerator}/{p_win_if_O1_corner.denominator}")
    
    term1 = p_O1_is_edge * p_win_if_O1_edge
    term2 = p_O1_is_corner * p_win_if_O1_corner
    print(f"P(Win) = {term1.numerator}/{term1.denominator} + {term2.numerator}/{term2.denominator}")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")

    # Return the answer in the specified format
    return f"<<<{total_win_prob.numerator}/{total_win_prob.denominator}>>>"

# Run the solver and print the final answer.
final_answer = solve_tic_tac_toe_probability()
print("\nFinal Answer:")
print(final_answer)
