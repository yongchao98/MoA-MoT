import math

# This script calculates the maximum distance from point A for an object shaped
# to create the strongest gravitational field at A, given a fixed volume of 1 m^3.

# The formula for this maximum distance 'd' is derived from calculus of variations
# and is given by d = (15 / (4 * pi))^(1/3).

# Define the numbers in the final equation
numerator = 15
denominator_coefficient = 4
power = 1/3

print("The formula for the furthest distance 'd' is:")
# Using an f-string to display the equation with its numerical components
print(f"d = ({numerator} / ({denominator_coefficient} * pi)) ^ ({power:.3f})")
print("-" * 30)

# Step 1: Calculate the value of the base of the exponent
base_value = numerator / (denominator_coefficient * math.pi)
print(f"Calculating the base value:")
print(f"{numerator} / ({denominator_coefficient} * {math.pi:.6f}) = {base_value:.6f}")
print("-" * 30)


# Step 2: Calculate the final result by applying the power (cube root)
final_distance = base_value ** power
print(f"Calculating the final distance (the cube root of the base):")
print(f"({base_value:.6f}) ^ ({power:.3f}) = {final_distance:.6f}")
print("-" * 30)

print(f"The furthest point on the surface of the playdough is {final_distance:.6f} meters from point A.")
