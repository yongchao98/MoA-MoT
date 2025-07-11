import fractions

def solve_yahtzee_problem():
    """
    Calculates and prints the solution to the Yahtzee probability problem.
    """
    # --- Strategy 1: Keep the '1' ---
    num_dice_to_roll_1 = 4
    target_face_1 = 1
    p1 = fractions.Fraction(1, 2**target_face_1)
    q1 = p1 + (1 - p1) * p1
    prob_1 = q1**num_dice_to_roll_1

    # --- Strategy 2: Keep the three '3's ---
    num_dice_to_roll_2 = 2
    target_face_2 = 3
    p3 = fractions.Fraction(1, 2**target_face_2)
    q3 = p3 + (1 - p3) * p3
    prob_2 = q3**num_dice_to_roll_2

    # --- Print results and the difference calculation ---
    print("This script calculates the difference in probabilities between two strategies in a game of Yahtzee.\n")
    print("--- Analysis of Each Strategy ---")
    print(f"Probability of success if keeping the '1':")
    print(f"P(Success | Strategy 1) = (P(roll '1' in 2 tries))^4 = ({q1})^4 = {prob_1}\n")
    
    print(f"Probability of success if keeping the three '3's:")
    print(f"P(Success | Strategy 2) = (P(roll '3' in 2 tries))^2 = ({q3})^2 = {prob_2}\n")

    # --- Calculate the difference ---
    difference = prob_1 - prob_2
    
    # Find common denominator for the final equation printout
    common_denom = difference.denominator
    p1_num = prob_1.numerator * (common_denom // prob_1.denominator)
    p2_num = prob_2.numerator * (common_denom // prob_2.denominator)
    
    print("--- Calculating the Difference ---")
    print("To find the difference, we subtract the probability of Strategy 2 from Strategy 1.")
    print("The final equation is:")
    print(f"{p1_num}/{common_denom} - {p2_num}/{common_denom} = {p1_num - p2_num}/{common_denom}")
    
    print(f"\nThe final difference as a fraction is {difference}.")
    print(f"As a decimal, the difference is approximately {float(difference):.4f}.")

solve_yahtzee_problem()
<<<1071/4096>>>