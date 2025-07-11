import numpy as np

# Plan:
# 1. Define the geometric constant 'a'.
# 2. Calculate the area moments of inertia I_ss and I_zz for the given cross-section.
#    The cross-section is a large square (3a x 3a) with two smaller squares (a x a) removed.
#    We use the parallel axis theorem: I = I_c + A*d^2.
#    I_total = I_large_square - I_hole_1 - I_hole_2.
# 3. Use I_ss and I_zz to calculate the parameters L and q0 based on the given formulas.
# 4. Use the derived symbolic formula F = (37 * q0 * L) / 270 to calculate the final force F.
# 5. Print the final calculation step and the result.

# Step 1: Define 'a'
a = 12**(1/4)
a_sq = np.sqrt(12)
a_fourth = 12

# Step 2: Calculate I_ss and I_zz
# Moment of inertia for a square of side 's' about its centroidal axis is s^4 / 12.

# For I_ss (about the s-axis, requires z-distances)
# Large square (3a x 3a), centered at origin
I_ss_large = (3*a) * (3*a)**3 / 12 = 81 * a_fourth / 12
# Hole 1 (a x a), centered at (s=a/2, z=-a)
I_ss_h1_c = a * a**3 / 12 = a_fourth / 12
A_h1 = a**2
d_z1 = -a
I_ss_h1 = I_ss_h1_c + A_h1 * d_z1**2 = a_fourth / 12 + a_sq * (-a)**2 = a_fourth / 12 + a_fourth
# Hole 2 (a x a), centered at (s=-a/2, z=a)
I_ss_h2_c = a_fourth / 12
A_h2 = a**2
d_z2 = a
I_ss_h2 = I_ss_h2_c + A_h2 * d_z2**2 = a_fourth / 12 + a_sq * a**2 = a_fourth / 12 + a_fourth

# Total I_ss
I_ss = I_ss_large - I_ss_h1 - I_ss_h2
# Simplified: I_ss = (81/12 - 2 * (1/12 + 1)) * a_fourth = (81/12 - 2 * 13/12) * a_fourth = (55/12) * a_fourth
I_ss_val = (55 / 12) * a_fourth

# For I_zz (about the z-axis, requires s-distances)
# Large square (3a x 3a)
I_zz_large = (3*a) * (3*a)**3 / 12 = 81 * a_fourth / 12
# Hole 1 (a x a), centered at (s=a/2, z=-a)
I_zz_h1_c = a**4 / 12 = a_fourth / 12
d_s1 = a/2
I_zz_h1 = I_zz_h1_c + A_h1 * d_s1**2 = a_fourth / 12 + a_sq * (a/2)**2 = a_fourth/12 + a_fourth/4
# Hole 2 (a x a), centered at (s=-a/2, z=a)
d_s2 = -a/2
I_zz_h2 = I_zz_h1_c + A_h2 * d_s2**2 = a_fourth/12 + a_sq * (-a/2)**2 = a_fourth/12 + a_fourth/4

# Total I_zz
I_zz = I_zz_large - I_zz_h1 - I_zz_h2
# Simplified: I_zz = (81/12 - 2 * (1/12 + 1/4)) * a_fourth = (81/12 - 2 * 4/12) * a_fourth = (73/12) * a_fourth
I_zz_val = (73 / 12) * a_fourth


# Step 3: Calculate L and q0
L = (30 * I_zz_val) / 73
q0 = (9 * I_ss_val) / 55

# Step 4: Calculate F
# F = (37 * q0 * L) / 270
F = (37 * q0 * L) / 270

# Step 5: Print the final answer
# The prompt requires printing the equation with the final numbers
print(f"I_ss = {I_ss_val}")
print(f"I_zz = {I_zz_val}")
print(f"L = (30 * {I_zz_val}) / 73 = {L}")
print(f"q0 = (9 * {I_ss_val}) / 55 = {q0}")
print(f"F = (37 * q0 * L) / 270")
print(f"F = (37 * {q0:.0f} * {L:.0f}) / 270 = {F:.2f}")
