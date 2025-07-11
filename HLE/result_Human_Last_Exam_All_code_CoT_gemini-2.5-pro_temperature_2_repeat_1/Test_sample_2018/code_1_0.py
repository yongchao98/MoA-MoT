# This script is designed to identify and present the value of a key
# computational factor from a specific scientific paper as described in the prompt.
# The factor is the "mushy region constant" from the original enthalpy-porosity
# method paper by Voller and Prakash (1987).

# The problem asks for the value originally used in the prior simulation-only work.
# Based on the paper, this value is 1.6 x 10^6.

# Define the components of the number in scientific notation.
base = 1.6
exponent = 6

# Calculate the final value of the computational factor.
computational_factor = base * (10 ** exponent)

# The prompt requires outputting the numbers in the final equation.
# This print statement shows the original components of the value.
print(f"Original computational factor equation: {base} * 10^{exponent}")

# Display the calculated result.
print(f"Resulting value: {computational_factor}")