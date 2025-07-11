from fractions import Fraction

# The probability of rolling side i is p(i) = 2^(-i)
def p(i):
    """Calculates the probability of rolling the number i."""
    return Fraction(1, 2**i)

# --- Strategy 1: Keep the '1' and aim for a Yahtzee of 1s. ---
# We have 4 dice to re-roll and need them to become '1'.
# We have two re-rolls.

# Probability of a die showing '1'
p_1 = p(1)

# The probability of a single die becoming '1' in at most two rolls is:
# P(1 on 1st try) + P(not 1 on 1st try) * P(1 on 2nd try)
prob_die_becomes_1 = p_1 + (1 - p_1) * p_1

# Since we need 4 dice to become '1', and the rolls are independent:
prob_S1 = prob_die_becomes_1 ** 4

# --- Strategy 2: Keep the three '3's and aim for a Yahtzee of 3s. ---
# We have 2 dice to re-roll and need them to become '3'.

# Probability of a die showing '3'
p_3 = p(3)

# The probability of a single die becoming '3' in at most two rolls is:
prob_die_becomes_3 = p_3 + (1 - p_3) * p_3

# Since we need 2 dice to become '3', and the rolls are independent:
prob_S2 = prob_die_becomes_3 ** 2

# --- Calculate the difference ---
difference = prob_S1 - prob_S2

# --- Print the explanation and result ---
print("This script calculates the difference in probabilities of achieving a Yahtzee between two strategies.")

print("\n--- Analysis of Strategy 1: Keep the '1' ---")
print("The goal is a Yahtzee of 1s. We need 4 dice to become '1's within two re-rolls.")
print(f"The probability of a single die becoming a '1' in one roll is p(1) = 1/2.")
print(f"The probability of one die becoming a '1' within two rolls is 1/2 + (1 - 1/2) * 1/2 = {prob_die_becomes_1.numerator}/{prob_die_becomes_1.denominator}.")
print(f"The total probability of success for Strategy 1 is ({prob_die_becomes_1.numerator}/{prob_die_becomes_1.denominator})^4 = {prob_S1.numerator}/{prob_S1.denominator}.")

print("\n--- Analysis of Strategy 2: Keep the three '3's ---")
print("The goal is a Yahtzee of 3s. We need 2 dice to become '3's within two re-rolls.")
print(f"The probability of a single die becoming a '3' in one roll is p(3) = 1/8.")
print(f"The probability of one die becoming a '3' within two rolls is 1/8 + (1 - 1/8) * 1/8 = {prob_die_becomes_3.numerator}/{prob_die_becomes_3.denominator}.")
print(f"The total probability of success for Strategy 2 is ({prob_die_becomes_3.numerator}/{prob_die_becomes_3.denominator})^2 = {prob_S2.numerator}/{prob_S2.denominator}.")

print("\n--- Difference in Probabilities ---")
print("The difference is P(Strategy 1) - P(Strategy 2).")
# To subtract fractions, find a common denominator. 4096 is a common denominator.
# P(Strategy 1) = 81/256 = (81*16)/(256*16) = 1296/4096
prob_S1_common_denom = prob_S1.limit_denominator(difference.denominator)
print(f"The final equation is: {prob_S1.numerator}/{prob_S1.denominator} - {prob_S2.numerator}/{prob_S2.denominator} = {prob_S1_common_denom.numerator}/{prob_S1_common_denom.denominator} - {prob_S2.numerator}/{prob_S2.denominator} = {difference.numerator}/{difference.denominator}")

print(f"\nThe final answer is {difference.numerator}/{difference.denominator}, which is approximately {float(difference):.4f}.")
