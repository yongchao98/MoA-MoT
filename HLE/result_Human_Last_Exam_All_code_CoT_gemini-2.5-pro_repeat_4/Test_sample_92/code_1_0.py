import fractions

# We need to solve the following system of linear equations for c1 and c2:
# 1) 9*c1 - c2 = 2
# 2) -c1 + 9*c2 = 6
#
# We are interested in c1, which is the probability of escape starting from bin 0.
#
# From equation (1), we can express c2 in terms of c1:
# c2 = 9*c1 - 2
#
# Now, we substitute this expression for c2 into equation (2):
# -c1 + 9*(9*c1 - 2) = 6
#
# Expanding this equation:
# -c1 + 81*c1 - 18 = 6
#
# This simplifies to the final equation for c1:
# 80*c1 = 24

c1_numerator = 24
c1_denominator = 80

# The final equation to solve for the probability (c1) is derived as follows.
# First, we state the system of equations derived from the problem's logic.
print("The system of equations for the probabilities c1 and c2 is:")
print("9*c1 - c2 = 2")
print("-1*c1 + 9*c2 = 6")
print("\nBy substituting the first equation into the second, we get the final equation for c1:")
# Outputting the numbers from the final equation
print(f"{c1_denominator} * c1 = {c1_numerator}")

# The solution for c1 is the probability we are looking for.
# We use the fractions module to get the result in its simplest form.
prob = fractions.Fraction(c1_numerator, c1_denominator)

print(f"\nThe probability that the marble escapes is c1 = {c1_numerator}/{c1_denominator}, which simplifies to {prob.numerator}/{prob.denominator}.")