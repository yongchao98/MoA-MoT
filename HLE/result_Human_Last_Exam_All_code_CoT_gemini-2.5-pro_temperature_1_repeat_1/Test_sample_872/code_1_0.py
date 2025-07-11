from fractions import Fraction

def solve_tic_tac_toe():
    """
    Calculates the maximum probability of winning Tic Tac Toe against a random opponent.
    The player goes first. A draw is considered a loss.
    """
    print("To maximize the chance of winning, the player (X) must choose the optimal first move.")
    print("We will analyze the case where X starts in the center, which is the optimal strategy.\n")

    # Step 1: X plays in the center. O has 8 squares left (4 corners, 4 edges).
    print("Step 1: Player X plays in the center.")
    print("The computer (O) has 8 empty squares to choose from randomly.")
    print("There are two cases for O's first move:")
    print("  A) O plays in a corner square (4 possibilities).")
    print("  B) O plays in an edge square (4 possibilities).\n")

    prob_O_corner = Fraction(4, 8)
    prob_O_edge = Fraction(4, 8)
    
    print(f"The probability O plays a corner is 4/8, which is {prob_O_corner}.")
    print(f"The probability O plays an edge is 4/8, which is {prob_O_edge}.\n")

    # --- Case A: O plays a corner ---
    print("--- Analyzing Case A: O plays a corner ---")
    print("X's best response creates an immediate threat. O has 6 available squares to respond.")
    
    # O's first chance to make a mistake
    prob_O_mistake1 = Fraction(5, 6)
    prob_O_correct1 = Fraction(1, 6)
    print(f"O must block X. The probability O fails to block (X wins) is {prob_O_mistake1.numerator}/{prob_O_mistake1.denominator}.")
    print(f"The probability O blocks correctly is {prob_O_correct1.numerator}/{prob_O_correct1.denominator}.")

    # If O blocks, X creates a second threat. O has 4 squares left.
    print("If O blocks, X creates a new threat. O has 4 available squares.")
    prob_O_mistake2 = Fraction(3, 4)
    prob_O_correct2 = Fraction(1, 4)
    print(f"  The probability O fails to block this second threat (X wins) is {prob_O_mistake2.numerator}/{prob_O_mistake2.denominator}.")
    print(f"  The probability O blocks this threat is {prob_O_correct2.numerator}/{prob_O_correct2.denominator}.")
    
    # If O blocks again, X creates a third threat. O has 2 squares left.
    print("If O blocks again, X creates a third threat. O has 2 available squares.")
    prob_O_mistake3 = Fraction(1, 2)
    prob_O_correct3 = Fraction(1, 2)
    print(f"    The probability O fails to block this third threat (X wins) is {prob_O_mistake3.numerator}/{prob_O_mistake3.denominator}.")
    print(f"    The probability O blocks, resulting in a draw (a loss for X), is {prob_O_correct3.numerator}/{prob_O_correct3.denominator}.")
    
    # Calculate total win probability for Case A
    win_if_O_correct3 = 0 # Draw is a loss
    win_after_2nd_block = prob_O_mistake3 * 1 + prob_O_correct3 * win_if_O_correct3
    win_after_1st_block = prob_O_mistake2 * 1 + prob_O_correct2 * win_after_2nd_block
    win_prob_case_A = prob_O_mistake1 * 1 + prob_O_correct1 * win_after_1st_block
    
    print("\nCalculating the total win probability for Case A:")
    print(f"P(Win | O corner) = P(O mistake 1) + P(O correct 1) * [P(O mistake 2) + P(O correct 2) * [P(O mistake 3) + P(O correct 3) * 0]]")
    print(f"P(Win | O corner) = {prob_O_mistake1} + {prob_O_correct1} * [{prob_O_mistake2} + {prob_O_correct2} * [{prob_O_mistake3}]]")
    print(f"P(Win | O corner) = {prob_O_mistake1} + {prob_O_correct1} * [{win_after_1st_block}] = {win_prob_case_A}\n")

    # --- Case B: O plays an edge ---
    print("--- Analyzing Case B: O plays an edge ---")
    print("X's best response again creates a threat that O must block from 6 available squares.")
    prob_O_mistake_B1 = Fraction(5, 6)
    prob_O_correct_B1 = Fraction(1, 6)
    print(f"The probability O fails to block (X wins) is {prob_O_mistake_B1.numerator}/{prob_O_mistake_B1.denominator}.")
    print(f"The probability O blocks correctly is {prob_O_correct_B1.numerator}/{prob_O_correct_B1.denominator}.")
    
    win_if_O_correct_B1 = 1 # Guaranteed win via a "fork"
    print("If O blocks correctly, X's next move creates a 'fork' (two ways to win), guaranteeing a win for X.")
    
    win_prob_case_B = prob_O_mistake_B1 * 1 + prob_O_correct_B1 * win_if_O_correct_B1
    
    print("\nCalculating the win probability for Case B:")
    print(f"P(Win|O edge) = P(O mistake) * 1 + P(O correct) * 1 (guaranteed win)")
    print(f"P(Win|O edge) = {prob_O_mistake_B1} + {prob_O_correct_B1} * 1 = {win_prob_case_B}\n")

    # --- Final calculation ---
    print("--- Calculating the Total Maximum Win Probability ---")
    total_win_prob = prob_O_corner * win_prob_case_A + prob_O_edge * win_prob_case_B
    part1 = prob_O_corner * win_prob_case_A
    part2 = prob_O_edge * win_prob_case_B
    
    print("Total P(Win) = P(O plays corner) * P(Win|O corner) + P(O plays edge) * P(Win|O edge)")
    print(f"Total P(Win) = ({prob_O_corner}) * ({win_prob_case_A}) + ({prob_O_edge}) * ({win_prob_case_B})")
    print(f"Total P(Win) = ({part1}) + ({part2})")
    print(f"Total P(Win) = {total_win_prob}")
    
    print("\nThis probability is higher than starting in a corner or on an edge.")
    print("Therefore, the maximum chance of winning is the value calculated above.")

    final_answer = total_win_prob
    print(f"\nThe final answer is {final_answer.numerator}/{final_answer.denominator}")

solve_tic_tac_toe()
<<<95/96>>>