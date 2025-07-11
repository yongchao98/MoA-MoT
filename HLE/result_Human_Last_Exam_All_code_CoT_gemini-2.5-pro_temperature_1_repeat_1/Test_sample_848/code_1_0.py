import math

# The problem asks for the integer part of 10^4 * C, where C is the limit
# C = lim_{N->inf} F(N)/ln(N).
# Our analysis shows that C = 2/ln(phi_13) + 2/ln(phi_5).
# The final equation to be calculated is:
# floor(10000 * (2/ln((13 + sqrt(165))/2) + 2/ln((5 + sqrt(21))/2)))

# Define the constants from the problem analysis
factor = 10000
m1 = 13
m2 = 5

# Calculate phi values, which are the growth rates of the two solution sequences.
# phi1 corresponds to m=13, and phi2 corresponds to m=5.
phi1 = (m1 + math.sqrt(m1**2 - 4)) / 2
phi2 = (m2 + math.sqrt(m2**2 - 4)) / 2

# Calculate the two terms that sum up to the limit C
term1 = 2 / math.log(phi1)
term2 = 2 / math.log(phi2)

# The limit C is the sum of these two terms
C = term1 + term2

# The value to be computed is 10^4 * C
final_value = factor * C

# The required answer is the integer part of this value
integer_part = int(final_value)

# Print the breakdown of the final equation as requested
print(f"The final equation to be calculated is:")
print(f"floor({factor} * (2/ln(({m1} + sqrt({m1**2 - 4}))/2) + 2/ln(({m2} + sqrt({m2**2 - 4}))/2)))")
print(f"\nBreaking down the calculation step-by-step:")
print(f"The first term in the sum corresponds to m = {m1}:")
print(f"  phi1 = ({m1} + sqrt({m1**2-4}))/2 = {phi1}")
print(f"  2 / ln(phi1) = {term1}")
print(f"\nThe second term in the sum corresponds to m = {m2}:")
print(f"  phi2 = ({m2} + sqrt({m2**2-4}))/2 = {phi2}")
print(f"  2 / ln(phi2) = {term2}")
print(f"\nThe limit C = (2 / ln(phi1)) + (2 / ln(phi2)) = {C}")
print(f"The expression inside the floor function = {factor} * {C} = {final_value}")
print(f"\nThe final integer part is: {integer_part}")