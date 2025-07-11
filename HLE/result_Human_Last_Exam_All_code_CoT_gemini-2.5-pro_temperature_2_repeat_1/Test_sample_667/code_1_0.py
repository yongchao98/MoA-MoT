from fractions import Fraction

# Step 1: Define the probability function for a single die roll
def p(i):
    """Calculates the probability of rolling the number i."""
    return Fraction(1, 2**i)

# Step 2: Calculate the probability of success for Strategy 1 (keeping the '1')
# The target number is 1. We re-roll 4 dice.
p_target_1 = p(1)
# Probability for a single die to become '1' in at most two re-rolls:
# P(success) = p_target_1 + (1 - p_target_1) * p_target_1
prob_one_die_becomes_1 = p_target_1 + (1 - p_target_1) * p_target_1
# Total probability for 4 dice to become '1's
prob_strategy_1 = prob_one_die_becomes_1 ** 4

# Step 3: Calculate the probability of success for Strategy 2 (keeping the '3's)
# The target number is 3. We re-roll 2 dice.
p_target_3 = p(3)
# Probability for a single die to become '3' in at most two re-rolls:
# P(success) = p_target_3 + (1 - p_target_3) * p_target_3
prob_one_die_becomes_3 = p_target_3 + (1 - p_target_3) * p_target_3
# Total probability for 2 dice to become '3's
prob_strategy_2 = prob_one_die_becomes_3 ** 2

# Step 4: Calculate the difference between the two probabilities
difference = prob_strategy_1 - prob_strategy_2

# Step 5: Print the final equation with all numbers
p1_num = prob_strategy_1.numerator
p1_den = prob_strategy_1.denominator
p2_num = prob_strategy_2.numerator
p2_den = prob_strategy_2.denominator
diff_num = difference.numerator
diff_den = difference.denominator

# To show the subtraction with a common denominator
common_den = p1_den * p2_den // difference.denominator
p1_adj_num = p1_num * (common_den // p1_den)
p2_adj_num = p2_num * (common_den // p2_den)


print(f"The probability of success by keeping the 1 is ({prob_one_die_becomes_1.numerator}/{prob_one_die_becomes_1.denominator})^4 = {p1_num}/{p1_den}.")
print(f"The probability of success by keeping the three 3s is ({prob_one_die_becomes_3.numerator}/{prob_one_die_becomes_3.denominator})^2 = {p2_num}/{p2_den}.")
print(f"The difference is {p1_num}/{p1_den} - {p2_num}/{p2_den} = {p1_adj_num}/{common_den} - {p2_adj_num}/{common_den} = {diff_num}/{diff_den}.")
<<<1071/4096>>>