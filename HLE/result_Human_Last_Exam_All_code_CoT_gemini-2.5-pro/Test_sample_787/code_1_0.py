import math

# Based on the number-theoretic analysis, we determine the exponents
# for the prime factors of the limit L.

# The exponent for the prime factor 2 is 10.
v2_L = 10

# The exponent for the prime factor 3 is 2.
v3_L = 2

# The exponent for the prime factor 5 is 1.
v5_L = 1

# The bases and their corresponding exponents
base1, exp1 = 2, v2_L
base2, exp2 = 3, v3_L
base3, exp3 = 5, v5_L

# Calculate the final limit
limit_value = (base1 ** exp1) * (base2 ** exp2) * (base3 ** exp3)

# Print the final equation with each number, as requested.
print(f"The limit of the sequence is calculated from the equation:")
print(f"{base1}^{exp1} * {base2}^{exp2} * {base3}^{exp3} = {limit_value}")
