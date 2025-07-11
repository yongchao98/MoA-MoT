# Based on the mathematical derivation, the exponents for the prime factors 2, 3, and 5
# of the limit L are determined to be 10, 2, and 1, respectively.

# Exponent for prime factor 2
e2 = 10
# Exponent for prime factor 3
e3 = 2
# Exponent for prime factor 5
e5 = 1

# Calculate the limit L
limit_g = (2**e2) * (3**e3) * (5**e5)

# Print the final equation and the result
print(f"The limit is given by the equation:")
print(f"2**{e2} * 3**{e3} * 5**{e5} = {limit_g}")
