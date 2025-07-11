import math

# Step 1: Define geometric and given parameters
a = 12**(1/4)
a_4 = 12.0 # a^4 = 12

# Step 2: Calculate the second moments of area I_zz and I_ss
# For the cross-section, which is a 3a x 3a square with two a x a squares removed.
# I_zz is about the vertical axis (s-axis), I_ss is about the horizontal axis (z-axis).
# Using the parallel axis theorem I = I_c + A*d^2
# Total I_zz = I_zz_main - 2 * I_zz_cutout
# I_zz = (81*a^4/12) - (a^4/3) - (a^4/3) = (73 * a^4) / 12
I_zz = (73 * a_4) / 12

# Total I_ss = I_ss_main - 2 * I_ss_cutout
# I_ss = (81*a^4/12) - (13*a^4/12) - (13*a^4/12) = (55 * a^4) / 12
I_ss = (55 * a_4) / 12

print(f"Calculated second moments of area:")
print(f"I_zz = {I_zz}")
print(f"I_ss = {I_ss}\n")


# Step 3: Calculate L and q0 based on the given data
# L = 30 * I_zz / 73
L = (30 * I_zz) / 73
# q0 = 9 * I_ss / 55
q0 = (9 * I_ss) / 55

print(f"Calculated beam parameters:")
print(f"L = {L}")
print(f"q0 = {q0}\n")

# Step 4: Use the derived formula to calculate F
# The relationship derived from setting the total deflection at x=3L/2 to zero is:
# F = (13 * q0 * L) / 270
F = (13 * q0 * L) / 270

# Step 5: Print the final equation with numerical values and the result
print("The final equation for the force F is:")
# Using format to show the numbers in the equation
print(f"F = (13 * {q0} * {L}) / 270")
print(f"F = {13 * q0 * L} / 270")
print(f"The required force F is: {F}")

# Final Answer Block
# The problem asks for the numerical value of the force F.
# So we wrap the final numerical answer in the required format.
# print(f'<<<{F}>>>')