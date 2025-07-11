from fractions import Fraction

# This function calculates the probability of a die becoming a target number
# within two rolls, given the single-roll probability 'p'.
def prob_success_two_rolls(p):
  """Calculates the probability of getting a target value within two rolls."""
  # This is equivalent to 1 - P(failure on roll 1 AND failure on roll 2)
  # P(failure) = 1 - p
  # P(failure on two rolls) = (1 - p) * (1 - p)
  return 1 - (1 - p)**2

# --- Strategy 1: Keep the '1' ---
# Goal: Get five '1's. We have one, need four more.
# The probability of rolling a '1' is p(1) = 2^(-1).
p_1 = Fraction(1, 2)
# The probability of a single die becoming '1' within two rolls.
prob_one_die_becomes_1 = prob_success_two_rolls(p_1)
# We need this to happen for all 4 dice we are re-rolling.
prob_strategy1 = prob_one_die_becomes_1 ** 4

# --- Strategy 2: Keep the three '3's ---
# Goal: Get five '3's. We have three, need two more.
# The probability of rolling a '3' is p(3) = 2^(-3).
p_3 = Fraction(1, 8)
# The probability of a single die becoming '3' within two rolls.
prob_one_die_becomes_3 = prob_success_two_rolls(p_3)
# We need this to happen for both 2 dice we are re-rolling.
prob_strategy2 = prob_one_die_becomes_3 ** 2

# --- Calculate the Difference ---
difference = prob_strategy1 - prob_strategy2

# --- Print the final equation ---
print("This problem asks for the difference in probabilities between two strategies.")
print("Strategy 1: Keep the '1' and re-roll four dice to get '1's.")
print(f"The probability of success for Strategy 1 is: {prob_strategy1.numerator}/{prob_strategy1.denominator}")
print("\nStrategy 2: Keep the three '3's and re-roll two dice to get '3's.")
print(f"The probability of success for Strategy 2 is: {prob_strategy2.numerator}/{prob_strategy2.denominator}")
print("\nThe difference is calculated as follows:")

print(f"\n{prob_strategy1.numerator}/{prob_strategy1.denominator} - {prob_strategy2.numerator}/{prob_strategy2.denominator} = {difference.numerator}/{difference.denominator}")