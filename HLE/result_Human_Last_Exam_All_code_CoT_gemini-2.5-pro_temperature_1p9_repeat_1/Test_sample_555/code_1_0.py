from fractions import Fraction

# The problem is to find the minimum expected number of rolls of a 7-sided die
# to generate a uniformly random digit from 0 to 9.

# 1. We always need at least two rolls. Let's call the baseline E_min = 2.
min_rolls = 2

# 2. After 2 rolls, there are 7*7 = 49 possible outcomes.
# We can map these 49 outcomes to our 10 target digits.
# A strategy is to take the integer value S from the two rolls (0-48) and see if
# the interval [S/49, (S+1)/49) lies within a [j/10, (j+1)/10) interval.
# The number of outcomes S for which we cannot decide is 9.
# These are S = {4, 9, 14, 19, 24, 29, 34, 39, 44}.
# The probability of being in one of these ambiguous states is 9/49.
prob_ambiguous_after_2_rolls = Fraction(9, 49)

# 3. If we are in an ambiguous state, we must roll again.
# It can be shown that for any ambiguous state, a single further roll
# will resolve the ambiguity with probability 6/7.
# The probability of needing yet another roll is 1/7.
prob_continue_from_ambiguous = Fraction(1, 7)

# 4. The expected number of *additional* rolls, given we are in an ambiguous state,
# follows a geometric distribution with p = (1 - 1/7) = 6/7.
# The expected value is 1/p.
expected_additional_rolls = 1 / (1 - prob_continue_from_ambiguous)

# 5. The total expected number of rolls is the baseline of 2, plus the
# probability of entering an ambiguous state multiplied by the expected
# additional rolls from that state.
# E = E_min + P(ambiguous) * E_additional
total_expected_rolls = min_rolls + prob_ambiguous_after_2_rolls * expected_additional_rolls

# Print out the components of the final calculation
print(f"The minimum number of rolls is always at least 2.")
print(f"The probability of needing more than 2 rolls is: {prob_ambiguous_after_2_rolls}")
print(f"The expected number of additional rolls given ambiguity is: {expected_additional_rolls}")
print("\nThe final calculation is:")
print(f"E = {min_rolls} + {prob_ambiguous_after_2_rolls.numerator}/{prob_ambiguous_after_2_rolls.denominator} * {expected_additional_rolls.numerator}/{expected_additional_rolls.denominator}")
print(f"E = {total_expected_rolls.numerator}/{total_expected_rolls.denominator}")
<<<31/14>>>