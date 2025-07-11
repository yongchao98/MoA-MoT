import math

# Define the constants from the problem statement and our derivation.
# The initial height H is given as 350g meters.
H_coefficient = 350
# The derived formula for time is T = sqrt(14 * H / g).
numerator_constant = 14

# We substitute H = 350g into the formula:
# T = sqrt((14 * (350 * g)) / g)
# The 'g' terms cancel out, leaving:
# T = sqrt(14 * 350)

# Calculate the value inside the square root
value_inside_sqrt = numerator_constant * H_coefficient

# Calculate the final time
time = math.sqrt(value_inside_sqrt)

# As requested, we print the numbers in the final equation step-by-step.
print(f"The final equation for time T is T = sqrt({numerator_constant} * {H_coefficient})")
print(f"T = sqrt({value_inside_sqrt})")
print(f"T = {time}")