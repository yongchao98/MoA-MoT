import math

# The length of segment AB is given.
ab_length = 18468

# Based on the geometric construction, the length of CP can be expressed
# as a multiple of the length of AB. The derived formula is:
# CP = AB * (sqrt(5) + 1) / 2
# The term (sqrt(5) + 1) / 2 is the golden ratio, often denoted by φ (phi).

# Let's perform the calculation step-by-step.
sqrt_5 = math.sqrt(5)
phi = (sqrt_5 + 1) / 2
cp_length = ab_length * phi

# Print out the equation with the substituted values
print(f"Let AB = {ab_length}")
print("The formula for the length of CP is: CP = AB * (sqrt(5) + 1) / 2\n")
print("Solving the equation:")
print(f"CP = {ab_length} * (math.sqrt(5) + 1) / 2")
print(f"CP = {ab_length} * ({sqrt_5} + 1) / 2")
print(f"CP = {ab_length} * {sqrt_5 + 1} / 2")
print(f"CP = {ab_length} * {phi}")
print(f"CP = {cp_length}\n")

# Print the final answer rounded to 4 decimal places.
print("The length of CP rounded to 4 decimal places is:")
print(f"{cp_length:.4f}")
print(f"<<<{cp_length:.4f}>>>")