from fractions import Fraction

def solve_tic_tac_toe_probability():
    """
    Calculates the maximum probability of winning Tic-Tac-Toe against a random opponent.
    The player goes first and a tie is considered a loss.
    """

    print("To maximize the chance of winning, your first move should be to take the center square.")
    print("The computer then chooses one of the remaining 8 squares at random.")
    print("Due to symmetry, we can analyze two main cases for the computer's first move:\n")

    # --- Case 1: Computer plays a corner square ---
    # There are 4 corner squares out of 8 total squares.
    p_o_corner = Fraction(4, 8)
    print(f"Case 1: Computer plays a corner square.")
    print(f"The probability of this is 4/8 = {p_o_corner}.")
    print("Your optimal response is to play on a non-adjacent edge to create an immediate threat.")
    print("For example, if you are 'X' at 5 and the computer is 'O' at 1, you play 'X' at 6.")
    print("You now have a line (5,6) and threaten to win at square 4.\n")

    # After your move, there are 6 empty squares.
    # You win if the computer does not play on your winning spot.
    p_win_on_t3_c = Fraction(5, 6)
    p_block_on_t3_c = Fraction(1, 6)
    print("The computer now has 6 empty squares to choose from.")
    print(f"You win immediately if the computer does not block your threat. Probability = {p_win_on_t3_c.numerator}/{p_win_on_t3_c.denominator}.")
    print(f"The game continues if the computer blocks. Probability = {p_block_on_t3_c.numerator}/{p_block_on_t3_c.denominator}.\n")
    
    # If the computer blocks, you must then block its threat and create a new one.
    # This leads to a situation where you can force another threat with 3 empty squares left.
    # The computer has a 1/3 chance of blocking, otherwise you win.
    p_win_on_t4_c = Fraction(2, 3)
    print("If the computer blocks, the game continues. Your next best move forces another threat.")
    print("The computer will have 3 squares left to choose from on its next turn.")
    print(f"You win if it fails to block this second threat. Probability = {p_win_on_t4_c.numerator}/{p_win_on_t4_c.denominator}.")
    print("If the computer blocks this second threat, the game results in a tie (which is a loss for you).\n")
    
    # Total probability of winning if the computer's first move is a corner.
    p_win_if_corner = p_win_on_t3_c + p_block_on_t3_c * p_win_on_t4_c
    print(f"So, the total probability of winning in Case 1 is:")
    print(f"P(win|corner) = P(win on turn 3) + P(block on T3) * P(win on turn 4)")
    print(f"P(win|corner) = {p_win_on_t3_c.numerator}/{p_win_on_t3_c.denominator} + ({p_block_on_t3_c.numerator}/{p_block_on_t3_c.denominator}) * ({p_win_on_t4_c.numerator}/{p_win_on_t4_c.denominator}) = {p_win_if_corner.numerator}/{p_win_if_corner.denominator}\n")

    # --- Case 2: Computer plays an edge square ---
    # There are 4 edge squares out of 8 total squares.
    p_o_edge = Fraction(4, 8)
    print("--------------------")
    print(f"Case 2: Computer plays an edge square.")
    print(f"The probability of this is 4/8 = {p_o_edge}.")
    print("Your optimal response is to play on an adjacent edge to create an immediate threat.")
    print("For example, if you are 'X' at 5 and the computer is 'O' at 2, you play 'X' at 4.")
    print("You now have a line (4,5) and threaten to win at square 6.\n")
    
    # The calculation for winning on turn 3 is identical.
    p_win_on_t3_e = Fraction(5, 6)
    p_block_on_t3_e = Fraction(1, 6)
    print("The computer now has 6 empty squares to choose from.")
    print(f"You win immediately if the computer does not block. Probability = {p_win_on_t3_e.numerator}/{p_win_on_t3_e.denominator}.")
    print(f"The game continues if the computer blocks. Probability = {p_block_on_t3_e.numerator}/{p_block_on_t3_e.denominator}.\n")

    # If the computer blocks, you again must block and create a new threat.
    # This leads to a situation where you can force another threat. The computer has 4 choices, 2 of which are losing.
    p_win_on_t4_e = Fraction(1, 2)
    print("If the computer blocks, your next best move forces another threat.")
    print(f"The computer has 4 choices, and a 1/2 chance of making a non-blocking move, allowing you to win.")
    print(f"Probability you win on the next turn = {p_win_on_t4_e.numerator}/{p_win_on_t4_e.denominator}.")
    print("If the computer blocks this second threat, the game results in a tie.\n")
    
    # Total probability of winning if the computer's first move is an edge.
    p_win_if_edge = p_win_on_t3_e + p_block_on_t3_e * p_win_on_t4_e
    print(f"So, the total probability of winning in Case 2 is:")
    print(f"P(win|edge) = {p_win_on_t3_e.numerator}/{p_win_on_t3_e.denominator} + ({p_block_on_t3_e.numerator}/{p_block_on_t3_e.denominator}) * ({p_win_on_t4_e.numerator}/{p_win_on_t4_e.denominator}) = {p_win_if_edge.numerator}/{p_win_if_edge.denominator}\n")
    
    # --- Final Calculation ---
    total_p_win = p_o_corner * p_win_if_corner + p_o_edge * p_win_if_edge
    print("--------------------")
    print("The total maximum probability of winning is the sum of the probabilities of these two cases:")
    print(f"Total P(Win) = P(O plays corner) * P(win|corner) + P(O plays edge) * P(win|edge)")
    print(f"Total P(Win) = ({p_o_corner.numerator}/{p_o_corner.denominator}) * ({p_win_if_corner.numerator}/{p_win_if_corner.denominator}) + ({p_o_edge.numerator}/{p_o_edge.denominator}) * ({p_win_if_edge.numerator}/{p_win_if_edge.denominator})")
    
    term1 = p_o_corner * p_win_if_corner
    term2 = p_o_edge * p_win_if_edge
    
    print(f"Total P(Win) = {term1.numerator}/{term1.denominator} + {term2.numerator}/{term2.denominator}")
    print(f"Total P(Win) = {total_p_win.numerator}/{total_p_win.denominator}")

    # Final Answer
    print(f"\nThus, the maximum chance of winning you can give yourself is {total_p_win.numerator}/{total_p_win.denominator}.")
    return total_p_win

final_answer = solve_tic_tac_toe_probability()
# <<<67/72>>>