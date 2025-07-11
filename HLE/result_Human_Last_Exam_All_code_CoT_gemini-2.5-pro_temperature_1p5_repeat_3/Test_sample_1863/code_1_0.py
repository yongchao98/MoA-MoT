import math

# Define the diameter of the quarter-sphere package in cm.
package_diameter = 250.0

# Calculate the theoretical maximum diameter of the inner sphere.
# The formula is derived from the geometry of the problem: d = D / (sqrt(3) + 1),
# where D is the package diameter and d is the inner sphere diameter.
sqrt_3_val = math.sqrt(3)
denominator = sqrt_3_val + 1
max_diameter_theoretical = package_diameter / denominator

# The diameters of available spheres are in increments of 0.01 cm.
# We must find the largest available diameter that is less than or equal to the
# theoretical maximum. This is achieved by rounding down (flooring) to two decimal places.
final_diameter = math.floor(max_diameter_theoretical * 100) / 100

print(f"A package in the shape of a quarter-sphere has a diameter of {package_diameter} cm.")
print("To find the maximum diameter (d) of a sphere that can fit inside, we use the formula:")
print("d = D / (sqrt(3) + 1)")
print("\n--- Calculation ---")
print(f"1. The numbers in the equation are:")
print(f"   D = {package_diameter}")
print(f"   sqrt(3) ≈ {sqrt_3_val}")
print("\n2. The equation with the numbers substituted is:")
print(f"   d = {package_diameter} / ({sqrt_3_val} + 1)")
print(f"   d = {package_diameter} / {denominator}")
print(f"   d ≈ {max_diameter_theoretical} cm")
print("\n3. Since the available spheres are in increments of 0.01 cm, we round down to the nearest available size.")
print(f"The maximum diameter of such a sphere is {final_diameter:.2f} cm.")