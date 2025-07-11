import math

# The value of r_0 is determined by the equation:
# r_0 = (226 + 49 * sqrt(2)) / 17

# Step 1: Calculate the value of the square root of 2
sqrt_2 = math.sqrt(2)
print(f"The value of sqrt(2) is: {sqrt_2}")

# Step 2: Calculate the numerator of the fraction
numerator = 226 + 49 * sqrt_2
print(f"The numerator (226 + 49 * sqrt(2)) is: {numerator}")

# Step 3: Define the denominator
denominator = 17
print(f"The denominator is: {denominator}")

# Step 4: Calculate the final value of r_0
r_0 = numerator / denominator
print(f"\nThe final equation is: r_0 = ({numerator}) / {denominator}")
print(f"The radial distance r_0 where the gravitational potential vanishes is: {r_0}")
