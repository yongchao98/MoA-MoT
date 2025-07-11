from fractions import Fraction

def solve_tic_tac_toe_prob():
    """
    Calculates the maximum probability of winning Tic Tac Toe against a random computer.

    The calculation is based on a game tree analysis from the optimal starting position.
    The player (X) goes first and aims to win (a tie is a loss).
    The computer (O) plays by choosing a random available square.

    The optimal strategy for the player is to start in the center.
    We then analyze the two symmetric cases for the computer's first move:
    1. Computer plays a corner.
    2. Computer plays an edge.
    """

    print("The player's optimal first move is to take the center square.")
    print("We analyze the two symmetric cases for the computer's first move.\n")

    # --- Case A: Computer plays a corner ---
    p_o_corner = Fraction(4, 8)

    # If O plays a corner, X must take the opposite corner. This threatens a win.
    # O has 6 squares left and must block one of two squares to survive.
    p_o_blocks_first_corner = Fraction(2, 6)
    p_o_fails_first_corner = Fraction(4, 6)

    # If O blocks the first threat, X must block O's subsequent threat.
    # This gives X a new threat that O must block from 4 remaining squares.
    p_o_blocks_second_corner = Fraction(1, 4)
    p_o_fails_second_corner = Fraction(3, 4)

    # If O blocks the second threat, the game continues. X can force a situation
    # where O has 2 moves left. One leads to a tie, the other to a win for X.
    # So, the win probability at this stage is 1/2.
    p_win_if_o_blocks_twice = Fraction(1, 2)

    # P(Win | O blocked first threat) = P(O fails second) * 1 + P(O blocks second) * P(Win|...)
    p_win_after_first_block = p_o_fails_second_corner + (p_o_blocks_second_corner * p_win_if_o_blocks_twice)

    # P(Win | O played a corner) = P(O fails first) * 1 + P(O blocks first) * P(Win|...)
    p_win_corner_case = p_o_fails_first_corner + (p_o_blocks_first_corner * p_win_after_first_block)

    # --- Case B: Computer plays an edge ---
    p_o_edge = Fraction(4, 8)

    # If O plays an edge, X's optimal move is a specific corner, creating a threat.
    # O has 6 squares left and must block.
    p_o_blocks_first_edge = Fraction(1, 6)
    p_o_fails_first_edge = Fraction(5, 6)

    # If O blocks, X can create a "fork" (two ways to win), guaranteeing a win.
    p_win_after_first_block_edge = Fraction(1, 1)

    # P(Win | O played an edge) = P(O fails first) * 1 + P(O blocks first) * 1
    p_win_edge_case = p_o_fails_first_edge + (p_o_blocks_first_edge * p_win_after_first_block_edge)


    # --- Final Calculation ---
    total_win_prob = (p_o_corner * p_win_corner_case) + (p_o_edge * p_win_edge_case)

    print("The final win probability is calculated as:")
    print("P(Win) = P(O plays corner) * P(Win | O plays corner) + P(O plays edge) * P(Win | O plays edge)")
    print(f"P(Win) = ({p_o_corner.numerator}/{p_o_corner.denominator}) * ({p_win_corner_case.numerator}/{p_win_corner_case.denominator}) + ({p_o_edge.numerator}/{p_o_edge.denominator}) * ({p_win_edge_case.numerator}/{p_win_edge_case.denominator})")
    
    term1 = p_o_corner * p_win_corner_case
    term2 = p_o_edge * p_win_edge_case
    
    print(f"P(Win) = ({term1.numerator}/{term1.denominator}) + ({term2.numerator}/{term2.denominator})")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")
    
    # The final answer format
    print("\nMaximum chance of winning:")
    print(f"<<<{total_win_prob.numerator}/{total_win_prob.denominator}>>>")

solve_tic_tac_toe_prob()