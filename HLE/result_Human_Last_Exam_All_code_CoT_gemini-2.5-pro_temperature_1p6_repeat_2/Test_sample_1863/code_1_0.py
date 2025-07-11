import math

# Define the given diameter of the quarter-sphere package.
package_diameter = 250.0

# Calculate the necessary components for the equation.
sqrt_of_2 = math.sqrt(2)
factor = sqrt_of_2 - 1
max_diameter = package_diameter * factor

# The formula for the maximum inner diameter 'd' is d = D * (sqrt(2) - 1),
# where D is the package diameter.
# The following lines show the calculation steps with the actual numbers.

print("Equation to find the maximum diameter:")
print(f"Maximum Diameter = Package Diameter * (sqrt(2) - 1)")
print("\nSubstituting the values:")
print(f"Maximum Diameter = {package_diameter} * ({sqrt_of_2} - 1)")
print("\nAfter calculating the factor:")
print(f"Maximum Diameter = {package_diameter} * {factor}")
print("\nFinal Result:")
print(f"Maximum Diameter = {max_diameter}")
