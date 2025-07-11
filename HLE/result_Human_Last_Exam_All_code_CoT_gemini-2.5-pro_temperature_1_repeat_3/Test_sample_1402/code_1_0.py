# The probability of Theo winning a single game is p_T = 1/6.
# The probability of Theo not winning a single game is 1 - p_T = 5/6.
# We want to find the probability that Theo does not win in the first four games.
# This is (5/6)^4.

# The base of the exponentiation
base_numerator = 5
base_denominator = 6

# The exponent
exponent = 4

# Calculate the result
result_numerator = base_numerator ** exponent
result_denominator = base_denominator ** exponent

# Print the final equation with all its components
print(f"The probability of Theo winning for the first time after at least five games is given by the equation:")
print(f"P = ({base_numerator}/{base_denominator})^{exponent}")
print(f"Calculating the numerator: {base_numerator}^{exponent} = {result_numerator}")
print(f"Calculating the denominator: {base_denominator}^{exponent} = {result_denominator}")
print(f"The final probability is the fraction: {result_numerator}/{result_denominator}")
print(f"As a decimal, the probability is approximately: {result_numerator/result_denominator:.4f}")