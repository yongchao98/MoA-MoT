import math
from fractions import Fraction

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# --- Strategy 1: Calculate the probability of success if we keep the three 3s ---

# Probability of rolling a 3 is 2^(-3)
p3 = Fraction(1, 8)
dice_to_roll_s1 = 2

# We have two rolls left.
# Probability of success on the 2nd roll (rolling two 3s)
prob_s1_roll2_success = p3**dice_to_roll_s1

# Probability of getting one 3 on the 2nd roll
prob_s1_roll2_one_3 = combinations(dice_to_roll_s1, 1) * p3**1 * (1 - p3)**1
# If we get one 3, we roll one die on the 3rd roll and need a 3
prob_s1_roll3_success_after_one_3 = prob_s1_roll2_one_3 * p3

# Probability of getting zero 3s on the 2nd roll
prob_s1_roll2_zero_3s = combinations(dice_to_roll_s1, 0) * p3**0 * (1 - p3)**2
# If we get zero 3s, we roll two dice on the 3rd roll and need two 3s
prob_s1_roll3_success_after_zero_3s = prob_s1_roll2_zero_3s * (p3**2)

# Total probability for the "keep 3s" strategy
prob_keep_3s = prob_s1_roll2_success + prob_s1_roll3_success_after_one_3 + prob_s1_roll3_success_after_zero_3s

# --- Strategy 2: Calculate the probability of success if we keep the one 1 ---

# Probability of rolling a 1 is 2^(-1)
p1 = Fraction(1, 2)
dice_to_roll_s2 = 4

# We have two rolls left.
# Probability of success on the 2nd roll (rolling four 1s)
prob_s2_roll2_success = p1**dice_to_roll_s2

# Calculate the probability of success on the 3rd roll based on the 2nd roll's outcome
prob_s2_roll3_success = 0
# k is the number of 1s we get on the 2nd roll (out of 4 dice)
for k in range(dice_to_roll_s2):  # k = 0, 1, 2, 3
    # Probability of getting k ones when rolling 4 dice
    prob_k_ones_on_roll2 = combinations(dice_to_roll_s2, k) * (p1**k) * ((1 - p1)**(dice_to_roll_s2 - k))
    
    # After getting k ones, we have k+1 ones total.
    # We need to roll 5 - (k+1) = 4-k more dice and get all 1s.
    dice_for_roll3 = dice_to_roll_s2 - k
    prob_success_on_roll3_given_k = p1**dice_for_roll3
    
    prob_s2_roll3_success += prob_k_ones_on_roll2 * prob_success_on_roll3_given_k

# Total probability for the "keep 1s" strategy
prob_keep_1s = prob_s2_roll2_success + prob_s2_roll3_success

# --- Final Output ---
difference = prob_keep_1s - prob_keep_3s

# Find a common denominator for printing the final equation
common_denominator = max(prob_keep_1s.denominator, prob_keep_3s.denominator, difference.denominator)
p1s_common = Fraction(prob_keep_1s.numerator * (common_denominator // prob_keep_1s.denominator), common_denominator)
p3s_common = Fraction(prob_keep_3s.numerator * (common_denominator // prob_keep_3s.denominator), common_denominator)
diff_common = Fraction(difference.numerator * (common_denominator // difference.denominator), common_denominator)

print("This problem asks for the difference in probability of achieving a Yahtzee between two strategies.")
print("\nStrategy 1: Keep the three 3s and aim for a Yahtzee of 3s.")
print(f"The probability of success for this strategy is {p3s_common.numerator}/{p3s_common.denominator}.")

print("\nStrategy 2: Keep the 1 and aim for a Yahtzee of 1s.")
print(f"The probability of success for this strategy is {p1s_common.numerator}/{p1s_common.denominator}.")

print("\nThe difference between the chance of getting a Yahtzee if you keep the 1 and if you keep the three 3s is:")
print(f"{p1s_common.numerator}/{p1s_common.denominator} - {p3s_common.numerator}/{p3s_common.denominator} = {diff_common.numerator}/{diff_common.denominator}")
