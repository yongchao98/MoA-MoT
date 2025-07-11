import math

# Plan:
# 1. Define the given geometric constant a^4.
# 2. Calculate the area moments of inertia, I_zz and I_ss, using the parallel axis theorem.
#    I_total = I_large_square - 2 * I_cutout_square
#    I_cutout_square = I_centroid + Area * distance^2
# 3. Use I_zz and I_ss to calculate L and q0.
# 4. Use the derived formula F = (13 * q0 * L) / 270 to find the force F.
# 5. Print the intermediate values and the final equation with numbers.

# Step 1: Define geometric constant
# a = 12^(1/4), so a^4 = 12
a_fourth = 12.0

# Step 2: Calculate moments of inertia
# For the main 3a x 3a square, I_centroid = (3a)(3a)^3 / 12 = 81 * a^4 / 12 = 27 * a_fourth / 4
I_main = (27.0 * a_fourth) / 4.0

# For I_zz (about z-axis, s-distance matters)
# Cutout is a x a, with center at (s, z). s_dist = +/- a/2.
# I_zz_cutout_centroid = a * a^3 / 12 = a_fourth / 12
# Area = a^2 = sqrt(a_fourth) * sqrt(a_fourth) - wait, this is not correct. a^2 is sqrt(144) = 12. Oh, a^2 is just sqrt(a_fourth). a^2 = sqrt(12). Let's recheck parallel axis theorem for a^2 d^2.
# Parallel axis term: A * d_s^2 = a^2 * (a/2)^2 = a^4 / 4 = a_fourth / 4
# My derivation was correct, a^2 * (a/2)^2 = a^4/4. Let's stick to a_fourth.
I_zz_cutout = (a_fourth / 12.0) + (a_fourth / 4.0) # I_c + A*d^2
I_zz = I_main - 2 * I_zz_cutout
I_zz_analytical = (73.0 / 12.0) * a_fourth # This is (27/4 - 2 * (1/12 + 1/4)) * a^4 = (27/4 - 2 * 4/12) * a^4 = (27/4 - 8/12) * a^4 = (81/12 - 8/12)*a^4 = 73/12 * a^4

# For I_ss (about s-axis, z-distance matters)
# Cutout is a x a, with center at (s, z). z_dist = +/- a.
# I_ss_cutout_centroid = a * a^3 / 12 = a_fourth / 12
# Parallel axis term: A * d_z^2 = a^2 * a^2 = a^4 = a_fourth
I_ss_cutout = (a_fourth / 12.0) + a_fourth # I_c + A*d^2
I_ss = I_main - 2 * I_ss_cutout
I_ss_analytical = (55.0 / 12.0) * a_fourth # This is (27/4 - 2 * (1/12 + 1)) * a^4 = (27/4 - 2 * 13/12) * a^4 = (27/4 - 26/12) * a^4 = (81/12 - 26/12)*a^4 = 55/12 * a^4

print(f"I_zz = {I_zz_analytical}")
print(f"I_ss = {I_ss_analytical}")

# Step 3: Calculate L and q0
L = (30.0 * I_zz_analytical) / 73.0
q0 = (9.0 * I_ss_analytical) / 55.0

print(f"L = {L}")
print(f"q0 = {q0}")

# Step 4: Calculate F
F = (13.0 * q0 * L) / 270.0

# Step 5: Print final equation and result
print("\nThe final equation for the force F is derived from setting the total deflection to zero:")
print("y_total = y_q + y_F = 0")
print("-(13 * q0 * L^4) / (240 * EI) + (9 * F * L^3) / (8 * EI) = 0")
print("This simplifies to: F = (13 * q0 * L) / 270")
print("\nSubstituting the calculated values:")
print(f"F = (13 * {q0} * {L}) / 270 = {F}")

print(f"\nThe required force F is {F}.")
print(f"<<<{F}>>>")