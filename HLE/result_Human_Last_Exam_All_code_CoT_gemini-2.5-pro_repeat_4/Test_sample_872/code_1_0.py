import fractions

def solve_tic_tac_toe_probability():
    """
    Calculates and explains the maximum probability of winning at Tic Tac Toe
    against a random opponent, assuming a tie is a loss.
    """
    print("To maximize the chance of winning, you should start by placing your 'X' on an edge square.")
    print("Let's number the board 1-9, and assume you play X at square 2 (a top edge square).\n")
    print(" 1 | X | 3 ")
    print("---|---|---")
    print(" 4 | 5 | 6 ")
    print("---|---|---")
    print(" 7 | 8 | 9 \n")
    print("The computer (O) will randomly choose one of the 8 remaining squares. We can group its move into 5 distinct cases by symmetry:\n")

    # Case 1: O plays center
    p_case1 = fractions.Fraction(1, 8)
    win_p_case1 = fractions.Fraction(5, 6)
    contrib1 = p_case1 * win_p_case1
    print(f"Case 1: O plays the center square (5).")
    print(f"  - Probability of this case: {p_case1.numerator}/{p_case1.denominator}")
    print(f"  - Your best response is to play at 1, threatening an immediate win at 3.")
    print(f"  - The computer must block at 3. It will do so with probability 1/6.")
    print(f"  - If it fails to block (probability 5/6), you win.")
    print(f"  - If it blocks, the game continues to a draw.")
    print(f"  - Your win probability in this case is {win_p_case1.numerator}/{win_p_case1.denominator}.")
    print(f"  - Contribution to total win probability: {p_case1.numerator}/{p_case1.denominator} * {win_p_case1.numerator}/{win_p_case1.denominator} = {contrib1.numerator}/{contrib1.denominator}\n")

    # Case 2: O plays an adjacent corner (1 or 3)
    p_case2 = fractions.Fraction(2, 8)
    win_p_case2 = fractions.Fraction(5, 6)
    contrib2 = p_case2 * win_p_case2
    print(f"Case 2: O plays an adjacent corner square (1 or 3).")
    print(f"  - Probability of this case: {p_case2.numerator}/{p_case2.denominator}")
    print(f"  - Your best response is to play at 5, threatening a win at 8.")
    print(f"  - If the computer fails to block at 8 (probability 5/6), you win.")
    print(f"  - If it blocks, the game leads to a draw.")
    print(f"  - Your win probability in this case is {win_p_case2.numerator}/{win_p_case2.denominator}.")
    print(f"  - Contribution to total win probability: {p_case2.numerator}/{p_case2.denominator} * {win_p_case2.numerator}/{win_p_case2.denominator} = {contrib2.numerator}/{contrib2.denominator}\n")

    # Case 3: O plays a non-adjacent corner (7 or 9)
    p_case3 = fractions.Fraction(2, 8)
    win_p_case3 = fractions.Fraction(1, 1) # This is a forced win
    contrib3 = p_case3 * win_p_case3
    print(f"Case 3: O plays a non-adjacent corner square (7 or 9).")
    print(f"  - Probability of this case: {p_case3.numerator}/{p_case3.denominator}")
    print(f"  - Your best response is to play at 5. You can then force a win by setting up a fork.")
    print(f"  - Your win probability in this case is {win_p_case3.numerator}/{win_p_case3.denominator}.")
    print(f"  - Contribution to total win probability: {p_case3.numerator}/{p_case3.denominator} * {win_p_case3.numerator}/{win_p_case3.denominator} = {contrib3.numerator}/{contrib3.denominator}\n")

    # Case 4: O plays an adjacent edge (4 or 6)
    p_case4 = fractions.Fraction(2, 8)
    win_p_case4 = fractions.Fraction(1, 1) # This is a forced win
    contrib4 = p_case4 * win_p_case4
    print(f"Case 4: O plays an adjacent edge square (4 or 6).")
    print(f"  - Probability of this case: {p_case4.numerator}/{p_case4.denominator}")
    print(f"  - Your best response is to play at 1. You can then force a win by setting up a fork.")
    print(f"  - Your win probability in this case is {win_p_case4.numerator}/{win_p_case4.denominator}.")
    print(f"  - Contribution to total win probability: {p_case4.numerator}/{p_case4.denominator} * {win_p_case4.numerator}/{win_p_case4.denominator} = {contrib4.numerator}/{contrib4.denominator}\n")
    
    # Case 5: O plays the opposite edge (8)
    p_case5 = fractions.Fraction(1, 8)
    win_p_case5 = fractions.Fraction(1, 1) # This is a forced win
    contrib5 = p_case5 * win_p_case5
    print(f"Case 5: O plays the opposite edge square (8).")
    print(f"  - Probability of this case: {p_case5.numerator}/{p_case5.denominator}")
    print(f"  - Your best response is to play at 1. You can then force a win by setting up a fork.")
    print(f"  - Your win probability in this case is {win_p_case5.numerator}/{win_p_case5.denominator}.")
    print(f"  - Contribution to total win probability: {p_case5.numerator}/{p_case5.denominator} * {win_p_case5.numerator}/{win_p_case5.denominator} = {contrib5.numerator}/{contrib5.denominator}\n")

    # Total probability
    total_prob = contrib1 + contrib2 + contrib3 + contrib4 + contrib5
    common_denominator = 48
    num1 = int(contrib1.numerator * (common_denominator / contrib1.denominator))
    num2 = int(contrib2.numerator * (common_denominator / contrib2.denominator))
    num3 = int(contrib3.numerator * (common_denominator / contrib3.denominator))
    num4 = int(contrib4.numerator * (common_denominator / contrib4.denominator))
    num5 = int(contrib5.numerator * (common_denominator / contrib5.denominator))
    
    print("The maximum chance of winning is the sum of the probabilities from all cases:")
    print(f"P(win) = {contrib1.numerator}/{contrib1.denominator} + {contrib2.numerator}/{contrib2.denominator} + {contrib3.numerator}/{contrib3.denominator} + {contrib4.numerator}/{contrib4.denominator} + {contrib5.numerator}/{contrib5.denominator}")
    print(f"       = {num1}/{common_denominator} + {num2}/{common_denominator} + {num3}/{common_denominator} + {num4}/{common_denominator} + {num5}/{common_denominator}")
    
    final_numerator = num1 + num2 + num3 + num4 + num5
    print(f"       = {final_numerator}/{common_denominator}")
    
    final_fraction = fractions.Fraction(final_numerator, common_denominator)
    print(f"\nThis fraction reduces to {final_fraction.numerator}/{final_fraction.denominator}.")

solve_tic_tac_toe_probability()
<<<15/16>>>