import math

# Step 1: Define the crystal properties and peak indices.
# Pseudocubic lattice parameter for NaMgH3 in Angstroms.
a = 3.91
# Miller indices (h, k, l) for the second major diffraction peak in a perovskite.
h, k, l = 2, 0, 0

# Step 2: Calculate the terms of the Q-space formula.
pi_val = math.pi
hkl_sq_sum = h**2 + k**2 + l**2
sqrt_hkl_sq_sum = math.sqrt(hkl_sq_sum)

# Step 3: Calculate the final Q value.
q_value = (2 * pi_val / a) * sqrt_hkl_sq_sum

# Step 4: Print the full calculation step-by-step.
print("The formula for the Q-space position of a peak in a cubic system is:")
print("Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)\n")
print("For the second major peak (2,0,0) of NaMgH3 (a = 3.91 Å):")

# The final code needs to output each number in the final equation.
print("Q = (2 * {:.5f} / {:.2f}) * sqrt({}^2 + {}^2 + {}^2)".format(pi_val, a, h, k, l))
print("Q = (2 * {:.5f} / {:.2f}) * sqrt({})".format(pi_val, a, hkl_sq_sum))
print("Q = (2 * {:.5f} / {:.2f}) * {:.5f}".format(pi_val, a, sqrt_hkl_sq_sum))
print("Q = {:.5f} * {:.5f}".format((2 * pi_val / a), sqrt_hkl_sq_sum))
print("Q = {:.4f} 1/Å".format(q_value))