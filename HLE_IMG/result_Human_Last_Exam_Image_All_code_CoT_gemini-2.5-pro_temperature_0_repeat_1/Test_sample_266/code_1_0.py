import math

# The side length of the red regular hexagon is given.
s = 3

# The formula for the area of a regular hexagon is (3 * sqrt(3) / 2) * s^2.
# Let's calculate the area.
area = (3 * math.sqrt(3) / 2) * (s**2)

# Print the steps of the calculation as requested.
print("Step 1: Identify the formula for the area of a regular hexagon.")
print("Formula: Area = (3 * sqrt(3) / 2) * s^2")
print("\nStep 2: Use the given side length s.")
print(f"s = {s}")
print("\nStep 3: Substitute the value of s into the formula.")
# The equation with the numbers plugged in
# Note: 3^2 = 9, and (3 * 9) / 2 = 27 / 2 = 13.5
term1 = 3
term2 = 3
term3 = 2
term4 = s
print(f"Area = ({term1} * sqrt({term2}) / {term3}) * {term4}^2")
print(f"Area = 13.5 * sqrt(3)")
print("\nStep 4: Calculate the final numerical value.")
print(f"Area ≈ {area:.2f}")
