import math

# Define the diameter of the quarter-sphere package
package_diameter = 250.0

# Calculate the radius of the package
package_radius = package_diameter / 2.0

# Calculate the theoretical maximum diameter of the inscribed sphere using the derived formula.
# The formula is d = package_diameter / (sqrt(3) + 1)
sqrt_3 = math.sqrt(3)
denominator = sqrt_3 + 1
theoretical_max_diameter = package_diameter / denominator

# Since available spheres have diameters in 0.01 cm increments, we must find the
# largest multiple of 0.01 that is less than or equal to the theoretical maximum.
# This is done by flooring the value to two decimal places.
final_diameter = math.floor(theoretical_max_diameter * 100) / 100

# --- Output the reasoning and the final equation with its numbers ---

print("Step 1: The radius of the quarter-sphere package (R) is the package diameter divided by 2.")
print(f"R = {package_diameter} / 2 = {package_radius} cm")
print("\nStep 2: The relationship between the package radius (R) and the inscribed sphere's radius (r) is r * (sqrt(3) + 1) = R.")
print("Solving for the diameter (d = 2r) in terms of the package diameter (D = 2R) gives:")
print("d = D / (sqrt(3) + 1)")
print("\nStep 3: Substitute the known values into the equation.")
print(f"d = {package_diameter} / (√3 + 1)")
print(f"d = {package_diameter} / ({sqrt_3} + 1)")
print(f"d = {package_diameter} / {denominator}")
print(f"The theoretical maximum diameter is {theoretical_max_diameter} cm.")

print("\nStep 4: Find the largest available diameter, which must be a multiple of 0.01 cm.")
print(f"The maximum available diameter is {final_diameter} cm.")

print(f"\n<<<{final_diameter}>>>")