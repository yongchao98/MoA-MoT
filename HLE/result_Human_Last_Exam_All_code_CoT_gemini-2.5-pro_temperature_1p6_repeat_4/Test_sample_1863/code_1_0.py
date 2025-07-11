import math

# The diameter of the quarter-sphere package.
package_diameter = 250.0

# The value of the square root of 3.
sqrt_3 = math.sqrt(3)

# The constant 1 in the denominator.
one = 1.0

# Calculate the theoretical maximum diameter using the formula: d = D / (sqrt(3) + 1)
theoretical_diameter = package_diameter / (sqrt_3 + one)

# Since available spheres have diameters in 0.01 cm increments, we must take the
# largest available size that is less than or equal to the theoretical maximum.
# This is equivalent to truncating the result to two decimal places.
final_diameter = math.floor(theoretical_diameter * 100) / 100.0

# Print the final equation with all its numerical components, as requested.
print(f"The maximum diameter is calculated using the equation: D / (sqrt(3) + 1)")
print(f"Substituting the values, we get:")
print(f"{package_diameter} / ({sqrt_3} + {one}) = {final_diameter}")
