import math

# Diameter of the quarter-sphere package in cm
package_diameter = 250

# The value of sqrt(3)
sqrt_3 = math.sqrt(3)

# Calculate the theoretical maximum diameter of the sphere that can fit inside.
# The formula is d = D / (sqrt(3) + 1), where d is the inner sphere's diameter
# and D is the package's diameter.
denominator = sqrt_3 + 1
max_diameter = package_diameter / denominator

# The available spheres have diameters in increments of 0.01 cm.
# We need to find the largest available diameter that is not greater than the calculated max_diameter.
# This can be done by truncating the result to two decimal places.
# We multiply by 100, take the integer part, and then divide by 100.
final_diameter = math.floor(max_diameter * 100) / 100

# Output the calculation steps and the final answer.
print("To find the maximum diameter of the sphere, we use the formula:")
print(f"d = D / (sqrt(3) + 1)")
print("Plugging in the values:")
print(f"d = {package_diameter} / ({sqrt_3:.4f} + 1)")
print(f"d = {package_diameter} / {denominator:.4f}")
print(f"d ≈ {max_diameter:.4f} cm")
print("\nSince the spheres come in diameter increments of 0.01 cm, the maximum diameter is the value above truncated to two decimal places.")
print(f"Maximum available diameter: {final_diameter:.2f} cm")

print(f"\n<<<The final answer is {final_diameter}>>>")