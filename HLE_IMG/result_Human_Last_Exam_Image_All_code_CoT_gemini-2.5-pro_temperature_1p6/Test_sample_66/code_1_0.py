import math

# Define the given values
R = 34  # Radius of the circle
L = 17  # Length of the tangent segment

# Step 1: Calculate the area of one of the large right-angled triangles
area_triangle = 0.5 * R * L

# Step 2: Calculate the area of the circular sector to be subtracted
# The angle alpha (in radians) is arctan(L/R)
angle_rad = math.atan(L / R)
# The area of the sector is (1/2) * R^2 * angle
area_sector = 0.5 * R**2 * angle_rad

# Step 3: Calculate the area of one shaded region
area_one_region = area_triangle - area_sector

# Step 4: Calculate the total shaded area (sum of the two regions)
total_area = 2 * area_one_region

# Print the formula and the final result
# The overall formula is 2 * (Area_triangle - Area_sector)
# This expands to 2 * ( (1/2)*R*L - (1/2)*R^2*arctan(L/R) ) which simplifies to R*L - R^2*arctan(L/R)
# Or 2 * Area_one_region
print("Step-by-step calculation:")
print(f"The radius (R) is {R}.")
print(f"The tangent length (L) is {L}.")
print(f"First, calculate the area of the large right-angled triangle: (1/2) * R * L = (1/2) * {R} * {L} = {area_triangle}")
print(f"Next, calculate the area of the circular sector to subtract. The angle is arctan(L/R) = arctan({L}/{R}) = {angle_rad:.4f} radians.")
print(f"The sector area is (1/2) * R^2 * angle = (1/2) * {R}^2 * {angle_rad:.4f} = {area_sector:.4f}")
print(f"The area of one shaded region is {area_triangle} - {area_sector:.4f} = {area_one_region:.4f}")
print("\nThe total area of the two shaded regions is double this amount.")
print(f"Total Area = 2 * ({area_triangle} - {area_sector:.4f})")
print(f"Total Area = {2*area_triangle} - 2 * {area_sector:.4f}")

# The problem asks to output the numbers in the final equation
# Final equation: Total Area = 2 * ( (1/2) * 34 * 17 - (1/2) * 34^2 * arctan(17/34) )
# Simplified: Total Area = 34 * 17 - 34^2 * arctan(17/34)
# = 578 - 1156 * arctan(0.5)
print("\nFinal Equation:")
print(f"Total Area = 2 * ( (1/2) * {R} * {L} - (1/2) * {R}**2 * arctan({L}/{R}) )")
print(f"Total Area = {R} * {L} - {R**2} * arctan({L/R})")
print(f"Total Area = {R*L} - {R**2} * {angle_rad:.4f}")
print(f"Total Area = {total_area:.4f}")
