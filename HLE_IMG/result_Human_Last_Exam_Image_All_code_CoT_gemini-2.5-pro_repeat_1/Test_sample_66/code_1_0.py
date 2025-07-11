import math

# Given values from the problem
R = 34  # Radius of the circular sector
L = 17  # Length of the tangent segment

# Area of one of the large right-angled triangles
area_triangle = 0.5 * R * L

# The angle of the sector in radians, derived from tan(alpha) = L/R
alpha = math.atan(L / R)

# Area of the circular sector to be subtracted from one triangle
area_sector = 0.5 * (R**2) * alpha

# The area of one shaded region (e.g., the blue part)
area_one_region = area_triangle - area_sector

# The total area is the sum of the two identical shaded regions
total_area = 2 * area_one_region

# The final equation is Total Area = 2 * (Triangle Area - Sector Area)
# Total Area = 2 * (289 - 578 * arctan(0.5))
# Total Area = 578 - 1156 * arctan(0.5)

# We print the components of the final equation and the result
term1 = 578
term2 = 1156
angle_val = math.atan(0.5)

print(f"The total area is calculated by the equation: {term1} - {term2} * arctan(0.5)")
print(f"The value of arctan(0.5) is approximately {angle_val:.4f} radians.")
print(f"The equation evaluates to: {term1} - {term2} * {angle_val:.4f} = {total_area:.4f}")
print("\nFinal Answer:")
print(f"{total_area}")
