import fractions

def solve_tic_tac_toe():
    """
    Calculates the maximum chance of winning Tic Tac Toe against a random opponent.

    The player (X) goes first and plays optimally to win. A tie is counted as a loss.
    The opponent (O) chooses a random unfilled square on each turn.

    The script calculates the win probability for the optimal opening move (a corner)
    and prints the detailed breakdown of the calculation.
    """
    print("To find the maximum chance of winning, we analyze the player's best opening move.")
    print("The optimal first move is to take a corner square. Let's calculate the win probability.\n")

    # Let X be the player and O be the random computer.
    # Turn 1: X takes a corner. 8 squares remain.

    # Turn 2: O plays randomly. We analyze two main cases based on O's move.
    prob_o_plays_center = fractions.Fraction(1, 8)
    prob_o_not_center = fractions.Fraction(7, 8)

    print(f"After X plays in a corner, O has 8 random choices.")
    print(f"Case A: O plays in the center. The probability of this is 1/8.")
    print(f"Case B: O does not play in the center. The probability of this is 7/8.\n")

    # --- Analysis of Case A: O plays center ---
    # Turn 3: X's optimal move is to play the opposite corner. 6 squares remain.
    # Turn 4: O has 6 random choices. These can be grouped into 2 types.
    prob_o_t4_corner = fractions.Fraction(2, 6)
    prob_o_t4_edge = fractions.Fraction(4, 6)

    # Subcase A1: O plays a corner on Turn 4.
    # X can block O's potential line, which also creates a fork for X. This is a guaranteed win.
    win_prob_A1 = fractions.Fraction(1, 1)

    # Subcase A2: O plays an edge on Turn 4.
    # For X to not win, O must make two consecutive perfect blocking moves.
    # The probability of O making the first correct block (on Turn 6) is 1/3.
    # The probability of O making the second correct block (on Turn 8) is 1/2.
    # The probability of O forcing a draw is (1/3) * (1/2) = 1/6.
    # Therefore, the probability of X winning in this subcase is 1 - 1/6 = 5/6.
    win_prob_A2 = fractions.Fraction(5, 6)
    
    # Total win probability for Case A
    win_prob_A = prob_o_t4_corner * win_prob_A1 + prob_o_t4_edge * win_prob_A2

    print("--- Analysis of Case A (O plays center) ---")
    print("X's best response is the opposite corner. Now O has 6 choices.")
    print(f" - If O plays one of the 2 other corners (Prob {prob_o_t4_corner}), X can force a win. Win prob = 1.")
    print(f" - If O plays one of the 4 edges (Prob {prob_o_t4_edge}), X wins unless O plays perfectly to force a draw.")
    print(f"   The probability of O forcing a draw is 1/6, so X's win probability here is {win_prob_A2}.")
    print(f"Total win probability for Case A = ({prob_o_t4_corner}) * {win_prob_A1} + ({prob_o_t4_edge}) * ({win_prob_A2}) = {win_prob_A}.\n")

    # --- Analysis of Case B: O does not play center ---
    # Turn 3: X takes the center. This creates an immediate threat.
    # Turn 4: O has 6 choices but must block X's threat to avoid an immediate loss.
    prob_o_fails_to_block = fractions.Fraction(5, 6)
    # If O fails to block, X wins immediately.
    # If O blocks, X can create a fork on the next move, guaranteeing a win.
    win_prob_B = fractions.Fraction(1, 1)
    
    print("--- Analysis of Case B (O does not play center) ---")
    print("X takes the center. X now threatens to win on the next turn.")
    print(f"O will fail to block this threat with probability {prob_o_fails_to_block}.")
    print("Even if O makes the 1 correct blocking move, X can force a win on the next turn.")
    print(f"So, the win probability for Case B is a guaranteed {win_prob_B}.\n")

    # --- Final Calculation ---
    total_win_prob = prob_o_plays_center * win_prob_A + prob_o_not_center * win_prob_B

    print("--- Final Calculation ---")
    print("The maximum win probability is found by combining these cases:")
    print("P(Win) = P(Win | Case A) * P(Case A) + P(Win | Case B) * P(Case B)")
    
    # Printing each number in the final equation
    print(f"P(Win) = ({win_prob_A.numerator}/{win_prob_A.denominator}) * ({prob_o_plays_center.numerator}/{prob_o_plays_center.denominator}) + ({win_prob_B.numerator}/{win_prob_B.denominator}) * ({prob_o_not_center.numerator}/{prob_o_not_center.denominator})")
    
    term1 = win_prob_A * prob_o_plays_center
    term2 = win_prob_B * prob_o_not_center
    
    print(f"P(Win) = ({term1.numerator}/{term1.denominator}) + ({term2.numerator}/{term2.denominator})")
    
    common_denominator = max(term1.denominator, term2.denominator)
    if common_denominator % term1.denominator != 0 or common_denominator % term2.denominator != 0:
        common_denominator = term1.denominator * term2.denominator

    term1_num = term1.numerator * (common_denominator // term1.denominator)
    term2_num = term2.numerator * (common_denominator // term2.denominator)

    print(f"P(Win) = {term1_num}/{common_denominator} + {term2_num}/{common_denominator}")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")

    # Return the final answer in the requested format
    global final_answer
    final_answer = f"{total_win_prob.numerator}/{total_win_prob.denominator}"

solve_tic_tac_toe()
print(f"\nThus, the maximum chance of winning is {final_answer}.")
<<<71/72>>>