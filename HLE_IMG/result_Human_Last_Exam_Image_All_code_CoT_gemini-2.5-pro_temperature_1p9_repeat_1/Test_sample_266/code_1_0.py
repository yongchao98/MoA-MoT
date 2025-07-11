import math

# Step 1: Define the side length of the red hexagon
s_red = 3

# Step 2: Calculate the area of the red hexagon
# Area_hexagon = (3 * sqrt(3) / 2) * s^2
area_red = (3 * math.sqrt(3) / 2) * (s_red ** 2)
print(f"The side length of the red hexagon is {s_red}.")
print(f"The area of the red hexagon is (3 * sqrt(3) / 2) * {s_red}^2 = {area_red:.4f}")


# Step 3: Use the area ratio to find the area of the white polygon.
# For this specific tiling, Area_white = (3/4) * Area_red.
area_ratio = 3/4
area_white = area_ratio * area_red
print(f"\nThe area of the white polygon is {area_ratio} times the area of the red hexagon.")
print(f"Area_white = {area_ratio} * {area_red:.4f} = {area_white:.4f}")

# Final Answer rounded to two decimal places
print(f"\nThe surface area of the smallest unit (the white shape) is approximately {area_white:.2f}.")
