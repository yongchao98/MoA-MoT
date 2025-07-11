import math

# Step 1: Define the given values from the geometric construction.
# 'a' is the length of the tangent leg of the right triangle.
# 'h' is the length of the hypotenuse.
a = 17
h = 34

# Step 2: Calculate the radius (R) of the circular sector using the Pythagorean theorem.
# R^2 + a^2 = h^2
R_squared = h**2 - a**2
R = math.sqrt(R_squared)

# Step 3: Determine the total central angle of the sector.
# The angle 2*theta is found from sin(2*theta) = a / h = 17 / 34 = 0.5.
# So, 2*theta is 30 degrees. The total angle of the sector shown is 4*theta.
total_angle_deg = 2 * 30 
total_angle_rad = math.radians(total_angle_deg)

# Step 4: Calculate the area of the kite formed by the two outer radii and their tangents.
# The kite's area is twice the area of one of the right triangles (with legs R and a).
# Area_Kite = 2 * (1/2 * R * a) = R * a
area_kite = R * a

# Step 5: Calculate the area of the circular sector.
# Area_Sector = (1/2) * R^2 * total_angle_rad
area_sector = 0.5 * R_squared * total_angle_rad

# Step 6: The total shaded area is the difference between the kite area and the sector area.
total_area = area_kite - area_sector

# Step 7: Print the final equation with all the numbers and the result.
# The base number for our calculations can be seen as 17*17 = 289.
base_val = a**2

print("The total shaded area is the area of the kite minus the area of the sector.")
print(f"Formula: Area = Area_Kite - Area_Sector")
print(f"The equation with the calculated numbers is:")
# The area of the kite is 289 * sqrt(3)
# The area of the sector is (1/2) * 867 * (pi/3) = 289 * pi / 2
print(f"Area = {base_val} * sqrt(3) - ({base_val} * pi / 2)")
print(f"Substituting the values of sqrt(3) and pi:")
print(f"Area = {base_val} * {math.sqrt(3):.4f} - ({base_val} * {math.pi:.4f} / 2)")
print(f"Area = {area_kite:.4f} - {area_sector:.4f}")
print(f"Final Total Area = {total_area:.4f}")