from fractions import Fraction

def calculate_yahtzee_prob(target_num, dice_to_roll):
    """
    Calculates the probability of turning a set of dice into a target number
    within two rolls.

    Args:
        target_num (int): The number we want the dice to show.
        dice_to_roll (int): The number of dice we are rolling.

    Returns:
        Fraction: The probability of success.
    """
    # Probability of rolling the target number in a single attempt
    p_single_roll = Fraction(1, 2**target_num)

    # Probability of a single die becoming the target number within two rolls
    # P(success) = P(success on 1st) + P(fail on 1st) * P(success on 2nd)
    prob_one_die_success = p_single_roll + (1 - p_single_roll) * p_single_roll

    # The total probability is for all dice to succeed independently
    total_prob = prob_one_die_success ** dice_to_roll
    return total_prob

# --- Strategy 1: Keep the '1' ---
prob_keep_1 = calculate_yahtzee_prob(target_num=1, dice_to_roll=4)

# --- Strategy 2: Keep the three '3's ---
prob_keep_3s = calculate_yahtzee_prob(target_num=3, dice_to_roll=2)

# --- Calculate the difference ---
difference = prob_keep_1 - prob_keep_3s

# --- Print the results ---
print("This script calculates the difference in probabilities for two Yahtzee strategies.")
print("\nStrategy 1: Keep the '1' and re-roll four dice.")
print(f"The probability of success is {prob_keep_1.numerator}/{prob_keep_1.denominator}.")

print("\nStrategy 2: Keep the three '3's and re-roll two dice.")
print(f"The probability of success is {prob_keep_3s.numerator}/{prob_keep_3s.denominator}.")

# To show the subtraction with a common denominator
common_denominator = prob_keep_1.denominator * prob_keep_3s.denominator / difference.denominator
prob_1_common_denom = prob_keep_1.numerator * (common_denominator / prob_keep_1.denominator)
prob_3_common_denom = prob_keep_3s.numerator * (common_denominator / prob_keep_3s.denominator)

print("\nCalculating the difference:")
print(f"P(Keep '1') - P(Keep '3's) = {prob_keep_1.numerator}/{prob_keep_1.denominator} - {prob_keep_3s.numerator}/{prob_keep_3s.denominator}")
print(f"                             = {int(prob_1_common_denom)}/{int(common_denominator)} - {int(prob_3_common_denom)}/{int(common_denominator)}")
print(f"                             = {difference.numerator}/{difference.denominator}")
print("\nThe final difference is {}/{}.".format(difference.numerator, difference.denominator))

# The required final answer format
# print(f"\n<<<{difference.numerator}/{difference.denominator}>>>")
# The question doesn't specify float or fraction, fraction is more precise.
# Let's give the decimal value as well.
final_value = float(difference.numerator) / float(difference.denominator)
print(f"As a decimal, this is approximately {final_value:.6f}.")

<<<1071/4096>>>