import math

# Step 1: Define the constants from the derived formula.
# The formula for the furthest distance R is R = (15 * V / (4 * pi))^(1/3).
# We are given the volume V = 1 m^3.
volume = 1.0
numerator = 15.0
denominator_constant = 4.0
pi_value = math.pi

# Step 2: Calculate the value of R.
# This R represents the maximum distance from point A to the surface of the playdough.
result = (numerator * volume / (denominator_constant * pi_value))**(1.0/3.0)

# Step 3: Print the final equation with the numbers plugged in, and the final answer.
# The instruction is to output each number in the final equation.
print("The equation for the furthest distance (R) is:")
print(f"R = ({numerator} / ({denominator_constant} * {pi_value}))^(1/3)")
print("\nCalculating this value gives:")
print(f"R = {result} meters")