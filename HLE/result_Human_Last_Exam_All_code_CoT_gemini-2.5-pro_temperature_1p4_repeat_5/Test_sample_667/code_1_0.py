from fractions import Fraction

def solve_yahtzee_problem():
    """
    Calculates the probability difference between two Yahtzee strategies
    with an infinitely-sided die.
    """
    # --- Strategy 1: Keep the '1' ---
    # We need 4 more '1's from 4 dice, with 2 rolls available.
    target_1 = 1
    dice_to_roll_1 = 4

    # Probability of rolling a '1'
    p_roll_1 = Fraction(1, 2**target_1)

    # Probability for a single die to become a '1' within two rolls
    # P(success) = P(1st roll is '1') + P(1st roll is not '1') * P(2nd roll is '1')
    p_success_one_die_1 = p_roll_1 + (1 - p_roll_1) * p_roll_1

    # Total probability for Strategy 1: all 4 dice must succeed
    prob_1 = p_success_one_die_1 ** dice_to_roll_1

    # --- Strategy 2: Keep the three '3's ---
    # We need 2 more '3's from 2 dice, with 2 rolls available.
    target_2 = 3
    dice_to_roll_2 = 2

    # Probability of rolling a '3'
    p_roll_2 = Fraction(1, 2**target_2)

    # Probability for a single die to become a '3' within two rolls
    p_success_one_die_2 = p_roll_2 + (1 - p_roll_2) * p_roll_2

    # Total probability for Strategy 2: both 2 dice must succeed
    prob_2 = p_success_one_die_2 ** dice_to_roll_2

    # --- Calculate the difference ---
    difference = prob_1 - prob_2

    # --- Print the results and the final equation ---
    print("Let P1 be the probability of success when keeping the '1'.")
    print(f"P1 = (chance for one die to become '1' in two rolls) ^ 4 = ({p_success_one_die_1}) ^ 4 = {prob_1}")
    print()
    print("Let P2 be the probability of success when keeping the three '3's.")
    print(f"P2 = (chance for one die to become '3' in two rolls) ^ 2 = ({p_success_one_die_2}) ^ 2 = {prob_2}")
    print()
    print("The difference is P1 - P2.")
    
    # To show the subtraction with a common denominator
    common_denominator = prob_2.denominator
    prob_1_common_denom_num = prob_1.numerator * (common_denominator // prob_1.denominator)
    prob_1_common_denom = Fraction(prob_1_common_denom_num, common_denominator)
    
    print("The final equation is:")
    print(f"{prob_1.numerator}/{prob_1.denominator} - {prob_2.numerator}/{prob_2.denominator} = {prob_1_common_denom.numerator}/{prob_1_common_denom.denominator} - {prob_2.numerator}/{prob_2.denominator} = {difference.numerator}/{difference.denominator}")
    
    # Returning the final numerical answer as a fraction string for the "answer" block
    return f"{difference.numerator}/{difference.denominator}"


final_answer = solve_yahtzee_problem()
# The final answer in the specified format
# The fraction 1071/4096 is approximately 0.261474609375
# For the format requested, just providing the fraction seems most appropriate
# and was explicitly calculated.
# The calculation is 81/256 - 225/4096 = 1296/4096 - 225/4096 = 1071/4096
# The result is 1071/4096.

# Let's put the answer in the special format
# The question doesn't ask for a float, so the fraction is the exact answer.
# >>>1071/4096
# However, the format is usually a number or letter, e.g. <<<C>>>, <<<9.8>>>
# Let's calculate the decimal value
final_decimal = 1071 / 4096
# print(f"\nFinal Answer as decimal: {final_decimal}")
# I will output the fraction as requested by the inner thoughts.
# The problem is a math problem, so a precise fraction is better than a rounded decimal.
print(f"\n<<<{final_answer}>>>")
