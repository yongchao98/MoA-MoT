import math

# Define the diameter of the sphere from which the package is cut.
package_sphere_diameter = 250.0

# The formula for the maximum diameter 'd' of a sphere that can fit
# inside a quarter-sphere of diameter 'D' is: d = D / (sqrt(3) + 1).
# Let's calculate the theoretical maximum diameter.
theoretical_diameter = package_sphere_diameter / (math.sqrt(3) + 1)

# The diameters of available spheres increment in steps of 0.01 cm.
# We must find the largest diameter that is less than or equal to the
# theoretical maximum. This is achieved by truncating to two decimal places.
final_diameter = math.floor(theoretical_diameter * 100) / 100

# Print the final equation with the numbers substituted.
print("The maximum theoretical diameter 'd_max' is calculated as follows:")
print(f"d_max = {package_sphere_diameter} / (√3 + 1)")
print(f"d_max = {package_sphere_diameter} / ({math.sqrt(3)} + 1)")
print(f"d_max ≈ {theoretical_diameter} cm")

# Print the final result after considering the available sizes.
print("\nSince the available sphere diameters are in increments of 0.01 cm, the largest possible sphere has a diameter of:")
print(f"{final_diameter:.2f} cm")
