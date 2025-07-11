import math

# The problem is to find the expected maximum earthquake magnitude, where magnitudes
# follow a Pareto(2) distribution and the number of observed earthquakes
# follows a LogSeries(1/2) distribution.

# Based on the derivation outlined in the plan, the final formula for the
# expected maximum magnitude (E) is:
# E = pi / (2 * ln(2)) - 1

print("The formula for the expected maximum earthquake magnitude (E) is:")
print("E = pi / (2 * ln(2)) - 1\n")

# Define the mathematical constants from the math library
pi_val = math.pi
ln2_val = math.log(2)

# Print the final equation with the numerical values of the constants
print("The final equation with the numerical values is:")
print(f"E = {pi_val} / (2 * {ln2_val}) - 1")

# Calculate the final result
result = pi_val / (2 * ln2_val) - 1

# Print the final calculated value
print("\nThe calculated value for the expected maximum magnitude is:")
print(result)