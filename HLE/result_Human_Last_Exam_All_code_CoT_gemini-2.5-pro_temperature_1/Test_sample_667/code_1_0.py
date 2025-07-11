from fractions import Fraction

# This script calculates the difference in probabilities for achieving a Yahtzee
# between two different strategies, based on the problem description.

# --- Part 1: Calculations for Strategy 1 (Keep the '1') ---

# The probability of rolling a '1' on the infinitely-sided die.
p1 = Fraction(1, 2**1)

# A die that is re-rolled has two chances to land on '1'.
# The probability is P(success on 1st try) + P(fail on 1st try) * P(success on 2nd try).
prob_one_die_becomes_1 = p1 + (1 - p1) * p1

# We need to turn four dice into '1's. The probability for this is the single-die
# probability raised to the power of 4.
prob_yahtzee_of_1s = prob_one_die_becomes_1 ** 4


# --- Part 2: Calculations for Strategy 2 (Keep the '3's) ---

# The probability of rolling a '3'.
p3 = Fraction(1, 2**3)

# The probability for one of the two re-rolled dice to become a '3' in two attempts.
prob_one_die_becomes_3 = p3 + (1 - p3) * p3

# We need to turn two dice into '3's. The probability is the single-die
# probability squared.
prob_yahtzee_of_3s = prob_one_die_becomes_3 ** 2


# --- Part 3: Print the explanation and final difference ---

print("This problem compares two strategies for getting a Yahtzee (all five dice the same).")
print("The initial roll is {1, 2, 3, 3, 3}, and the probability of rolling side 'i' is 2**(-i).")
print("We have two re-rolls available for any dice we don't keep.\n")

print("--- Strategy 1: Keep the '1' and aim for a Yahtzee of 1s ---")
print("The probability for one die to become a '1' in two rolls is:")
print(f"P(die becomes 1) = P(1) + (1 - P(1)) * P(1) = {p1} + (1 - {p1}) * {p1} = {prob_one_die_becomes_1}")
print("Since we need to turn four dice into '1's, the total probability is:")
print(f"P(Yahtzee of 1s) = ({prob_one_die_becomes_1})^4 = {prob_yahtzee_of_1s}\n")

print("--- Strategy 2: Keep the three '3's and aim for a Yahtzee of 3s ---")
print("The probability for one die to become a '3' in two rolls is:")
print(f"P(die becomes 3) = P(3) + (1 - P(3)) * P(3) = {p3} + (1 - {p3}) * {p3} = {prob_one_die_becomes_3}")
print("Since we need to turn two dice into '3's, the total probability is:")
print(f"P(Yahtzee of 3s) = ({prob_one_die_becomes_3})^2 = {prob_yahtzee_of_3s}\n")

print("--- Final Difference ---")
print("The difference between the probabilities of the two strategies is:")

# To show the subtraction with a common denominator
difference = prob_yahtzee_of_1s - prob_yahtzee_of_3s
common_denominator = difference.denominator
p1_common_num = prob_yahtzee_of_1s.numerator * (common_denominator // prob_yahtzee_of_1s.denominator)
p3_common_num = prob_yahtzee_of_3s.numerator * (common_denominator // prob_yahtzee_of_3s.denominator)

# Print the final equation with all numbers
print(f"({prob_one_die_becomes_1})^4 - ({prob_one_die_becomes_3})^2 = {prob_yahtzee_of_1s} - {prob_yahtzee_of_3s} = {p1_common_num}/{common_denominator} - {p3_common_num}/{common_denominator} = {difference}")