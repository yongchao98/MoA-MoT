# This script calculates the limiting probability P(n) as n approaches infinity.
# Based on the analysis, the limiting value is interpreted as the average probability
# over its oscillation cycle, which simplifies to the case of 8 independent events.

# The probability of the collapse condition NOT being met for a single base is 1 - 1/8 = 7/8.
# There are 8 bases, and in the averaged limit, we can treat them as independent.
# The final probability is (7/8)^8.

# The numerator of the fraction is 7^8.
numerator = 7**8

# The denominator of the fraction is 8^8.
denominator = 8**8

print("The closed-form expression for the limiting probability is (7^8) / (8^8).")
print(f"The numerator of the final equation is 7^8 = {numerator}")
print(f"The denominator of the final equation is 8^8 = {denominator}")

# Calculate the final numerical value.
limiting_probability = numerator / denominator
print(f"The limiting value of P(n) is: {limiting_probability}")