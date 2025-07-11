import fractions

def solve_tic_tac_toe_probability():
    """
    Calculates and explains the maximum win probability in Tic-Tac-Toe against a random player.

    The optimal strategy is to start in the center. The calculation proceeds as follows:
    1. You (X) place your mark in the center square.
    2. The computer (O) has 8 empty squares to choose from randomly. These are categorized into
       4 corners and 4 edges.

    Case A: The computer plays a corner.
       - You respond by taking the opposite corner. This creates a situation where you can win if the
         computer fails to block on its next move.
       - There will be 6 squares left for the computer's second move. 4 of these moves will result
         in you being able to force a win. 2 of them will block you, leading to a draw (a loss).
       - So, your win probability in this case is 4/6.

    Case B: The computer plays an edge.
       - You respond by taking a specific corner that sets up a threat.
       - Even if the computer blocks your immediate threat (1 out of 6 moves), you can make a move
         that creates an unstoppable "fork" (two ways to win on the next turn).
       - Therefore, your win probability in this case is 1.

    We combine these cases to get the total maximum probability of winning.
    """

    # Probabilities for the computer's first move after you take the center.
    # There are 8 squares left: 4 corners, 4 edges.
    p_comp_corner = fractions.Fraction(4, 8)
    p_comp_edge = fractions.Fraction(4, 8)

    # Win probability if the computer plays a corner first.
    # After your response, 6 squares are left. 4 of the computer's possible moves
    # will lead to your win.
    p_win_if_comp_corner = fractions.Fraction(4, 6)

    # Win probability if the computer plays an edge first.
    # Your optimal response leads to a forced win.
    p_win_if_comp_edge = fractions.Fraction(1, 1)

    # Total win probability is the sum of probabilities of the two disjoint cases.
    total_p_win = (p_comp_corner * p_win_if_comp_corner) + (p_comp_edge * p_win_if_comp_edge)

    print("To achieve the maximum chance of winning, you must start by placing your 'X' in the center square.")
    print("The computer then chooses an unfilled square randomly. We analyze the two possible cases:\n")
    print("1. The computer picks a corner square.")
    print("2. The computer picks an edge square.\n")

    print("The calculation is as follows:")
    print("P(Win) = P(Computer plays Corner) * P(Win | Given Computer plays Corner) + P(Computer plays Edge) * P(Win | Given Computer plays Edge)\n")

    print(f"The probability the computer plays a corner is {p_comp_corner.numerator}/{p_comp_corner.denominator}.")
    print(f"In this case, your chance of winning is {p_win_if_comp_corner.numerator}/{p_win_if_comp_corner.denominator}.\n")

    print(f"The probability the computer plays an edge is {p_comp_edge.numerator}/{p_comp_edge.denominator}.")
    print(f"In this case, your chance of winning is {p_win_if_comp_edge.numerator}/{p_win_if_comp_edge.denominator}.\n")

    term1 = p_comp_corner * p_win_if_comp_corner
    term2 = p_comp_edge * p_win_if_comp_edge
    
    # We simplify the terms before printing the equation to make it clearer.
    p_comp_corner_s = fractions.Fraction(p_comp_corner).limit_denominator()
    p_win_if_comp_corner_s = fractions.Fraction(p_win_if_comp_corner).limit_denominator()
    p_comp_edge_s = fractions.Fraction(p_comp_edge).limit_denominator()

    print("So, the final equation is:")
    print(f"P(Win) = ({p_comp_corner_s.numerator}/{p_comp_corner_s.denominator}) * ({p_win_if_comp_corner_s.numerator}/{p_win_if_comp_corner_s.denominator}) + ({p_comp_edge_s.numerator}/{p_comp_edge_s.denominator}) * ({p_win_if_comp_edge.numerator}/{p_win_if_comp_edge.denominator})")
    print(f"P(Win) = ({term1.numerator}/{term1.denominator}) + ({term2.numerator}/{term2.denominator})")
    
    # To show the addition with a common denominator
    term1_final = fractions.Fraction(term1).limit_denominator()
    term2_final = fractions.Fraction(term2).limit_denominator()
    
    print(f"P(Win) = {term1_final.numerator}/{term1_final.denominator} + {term2_final.numerator}/{term2_final.denominator}")
    print(f"P(Win) = {total_p_win.numerator}/{total_p_win.denominator}\n")

    print(f"The maximum chance of winning is {total_p_win.numerator}/{total_p_win.denominator}.")
    print(f"<<<{total_p_win.numerator}/{total_p_win.denominator}>>>")


if __name__ == '__main__':
    solve_tic_tac_toe_probability()