from fractions import Fraction

def calculate_yahtzee_probability(target_number, dice_to_roll):
    """
    Calculates the probability of getting a specific Yahtzee given a number of dice to roll.
    
    Args:
        target_number (int): The number we are aiming for (e.g., 1 for five 1s).
        dice_to_roll (int): The number of dice we need to change to the target_number.
        
    Returns:
        A tuple containing:
        - prob_one_die (Fraction): The probability of one die becoming the target number.
        - total_prob (Fraction): The total probability of all dice becoming the target number.
    """
    # Probability of rolling the target number in a single roll
    p_k = Fraction(1, 2**target_number)
    
    # Probability of one die becoming the target number within two rolls
    # This is P(success on 1st roll) + P(fail on 1st AND success on 2nd)
    prob_one_die = p_k + (1 - p_k) * p_k
    
    # Total probability for all dice to become the target number
    total_prob = prob_one_die ** dice_to_roll
    
    return prob_one_die, total_prob

# --- Strategy 1: Keep the '1' ---
# We have one '1', so we need to roll the other 4 dice to get '1's.
target_1 = 1
dice_to_roll_1 = 4
prob_die_is_1, p1 = calculate_yahtzee_probability(target_1, dice_to_roll_1)

# --- Strategy 2: Keep the three '3's ---
# We have three '3's, so we need to roll the other 2 dice to get '3's.
target_3 = 3
dice_to_roll_3 = 2
prob_die_is_3, p2 = calculate_yahtzee_probability(target_3, dice_to_roll_3)

# --- Calculate the Difference ---
difference = p1 - p2

# --- Print the results and the final equation ---

print("Step 1: Calculate the probability of success for keeping the '1'.")
print(f"The probability of a die becoming a '1' in two rolls is {prob_die_is_1.numerator}/{prob_die_is_1.denominator}.")
print(f"The total probability for four dice is ({prob_die_is_1.numerator}/{prob_die_is_1.denominator})^{dice_to_roll_1}, which is {p1.numerator}/{p1.denominator}.\n")

print("Step 2: Calculate the probability of success for keeping the three '3's.")
print(f"The probability of a die becoming a '3' in two rolls is {prob_die_is_3.numerator}/{prob_die_is_3.denominator}.")
print(f"The total probability for two dice is ({prob_die_is_3.numerator}/{prob_die_is_3.denominator})^{dice_to_roll_3}, which is {p2.numerator}/{p2.denominator}.\n")

print("Step 3: Calculate the difference between the two probabilities.")
# To subtract, we find a common denominator.
common_denominator = difference.denominator
p1_common_num = p1.numerator * (common_denominator // p1.denominator)

print("The final equation for the difference is:")
print(f"({p1.numerator}/{p1.denominator}) - ({p2.numerator}/{p2.denominator}) = ({p1_common_num}/{common_denominator}) - ({p2.numerator}/{common_denominator}) = ({p1_common_num - p2.numerator})/{common_denominator}")

print(f"\nThe exact difference is {difference.numerator}/{difference.denominator}.")
print(f"As a decimal, this is approximately {float(difference):.4f}.")

# The final answer in the required format
final_answer = float(difference)
# <<<1071/4096>>> # I'll output the decimal representation as per example in prompt.
final_answer_decimal = round(final_answer, 4)
# print(f'<<<{final_answer_decimal}>>>') # Hidden from user to avoid confusion. Final format at end.