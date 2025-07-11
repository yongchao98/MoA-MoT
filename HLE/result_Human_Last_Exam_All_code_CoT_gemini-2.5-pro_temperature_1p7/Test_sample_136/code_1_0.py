import math

# This script calculates the limit of the conditional probability P[x_0 not visited | 0 not visited]
# The final formula for the limit is exp(-(ln(2)^2 / (4 * pi)))

# 1. Calculate the components of the formula
ln_2 = math.log(2)
ln_2_squared = ln_2 ** 2
four_pi = 4 * math.pi

# 2. Calculate the exponent
exponent = -ln_2_squared / four_pi

# 3. Calculate the final result
limit_probability = math.exp(exponent)

# 4. Print the final equation with all its numerical components
print("The limit is calculated using the formula: exp(-(ln(2)^2 / (4 * pi)))")
print("\nStep-by-step calculation:")
print(f"Value of ln(2): {ln_2}")
print(f"Value of (ln(2))^2: {ln_2_squared}")
print(f"Value of 4*pi: {four_pi}")
print(f"The exponent value (-(ln(2)^2 / (4 * pi))): {exponent}")
print(f"\nThe final limit of the conditional probability is: {limit_probability}")
