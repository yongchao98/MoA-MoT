import numpy as np

# Based on the step-by-step derivation, the parameters R and phi for the
# deferent-epicycle model are determined by the Fourier series representation
# of the object's motion.

# The analysis shows that the two dominant components of the motion correspond to
# the Fourier series terms with indices k=1 and k=-3.

# The radius of the deferent (the larger circle) is determined by the magnitude
# of the k=1 coefficient, |c_1|.
# The radius of the epicycle (the smaller circle) is determined by the magnitude
# of the k=-3 coefficient, |c_{-3}|.

# The frequencies are given by the indices k:
# Deferent frequency: Omega_def = 1 * omega_0
# Epicycle frequency: Omega_epi = -3 * omega_0
# (where omega_0 is the fundamental frequency of the orbit around the square)

# R is the ratio of the radii, R = |c_1| / |c_{-3}|.
# The analytic derivation results in the expression R = 18 / (3*pi - 2).
R_numerator = 18
R_denominator_coeff = 3
R_denominator_const = -2
R = R_numerator / (R_denominator_coeff * np.pi + R_denominator_const)

# phi is the ratio of the frequencies, phi = Omega_epi / Omega_def.
# phi = (-3 * omega_0) / (1 * omega_0) = -3.
phi_numerator = -3
phi_denominator = 1
phi = phi_numerator / phi_denominator

# The final code outputs the formula for each number and its calculated value.
print("To find the ordered pair (R, phi):")
print("\n1. We calculate R, the ratio of the deferent radius to the epicycle radius.")
print(f"The formula for R is: {R_numerator} / ({R_denominator_coeff} * pi + ({R_denominator_const}))")
print(f"   R = {R}")

print("\n2. We calculate phi, the ratio of the epicycle frequency to the deferent frequency.")
print(f"The formula for phi is: {phi_numerator} / {phi_denominator}")
print(f"   phi = {phi}")

print(f"\nThe resulting ordered pair (R, phi) is approximately ({R:.4f}, {phi:.4f}).")
