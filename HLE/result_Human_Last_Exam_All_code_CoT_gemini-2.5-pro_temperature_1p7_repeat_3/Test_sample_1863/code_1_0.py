import math

# The diameter of the quarter-sphere package in cm.
package_diameter = 250.0

# --- Calculation ---
# The formula for the maximum diameter (d) of a sphere that fits in a
# quarter-sphere of diameter (D) is d = D / (sqrt(3) + 1).

# 1. Calculate the value of the square root of 3.
sqrt_3 = math.sqrt(3)

# 2. Calculate the denominator for the formula.
denominator = sqrt_3 + 1

# 3. Calculate the theoretical maximum diameter.
theoretical_diameter = package_diameter / denominator

# 4. The available spheres have diameters in 0.01 cm increments.
#    We find the largest diameter that fits by truncating the theoretical
#    result to two decimal places.
final_diameter = math.floor(theoretical_diameter * 100) / 100

# --- Output ---
print("To find the maximum diameter of the sphere, we use the formula: d = D / (sqrt(3) + 1)")
print("Here, 'D' is the diameter of the quarter-sphere package and 'd' is the diameter of the sphere to fit inside.")
print("\nHere is the step-by-step calculation with the given numbers:")
print(f"The final equation is: d = {package_diameter} / (math.sqrt(3) + 1)")
print(f"First, we calculate the denominator: ( {sqrt_3} + 1 ) = {denominator}")
print(f"Next, we solve for 'd': d = {package_diameter} / {denominator}")
print(f"The theoretical maximum diameter 'd' is approximately {theoretical_diameter} cm.")
print("\nSince the available spheres have diameters in increments of 0.01 cm, the maximum diameter must be less than or equal to the theoretical value.")
print(f"By truncating to two decimal places, we find the maximum practical diameter.")
print(f"\nThe maximum diameter of a sphere that can fit in the package is {final_diameter} cm.")