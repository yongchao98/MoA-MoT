import math

# Step 1: Define the geometric parameter 'a'
a = 12**(1/4)
# We know a^4 = 12, which simplifies calculations
a4 = 12.0

# Step 2: Calculate the moments of inertia I_ss and I_zz
# I_zz = moment of inertia about the s-axis
# I_ss = moment of inertia about the z-axis
# I_main = I for the 3a x 3a square about its centroid
# I_cutout_c = I for the a x a square about its centroid
I_main = (3*a) * (3*a)**3 / 12
I_cutout_c = a * a**3 / 12

# For I_zz, the distance 'd' is in the z-direction. d1 = -a, d2 = a
# I_zz = I_main - 2 * (I_cutout_c + Area_cutout * d^2)
Area_cutout = a**2
d_z = a
I_zz = I_main - 2 * (I_cutout_c + Area_cutout * d_z**2)
I_zz = (81/12)*a4 - 2 * ( (1/12)*a4 + a*a*a*a )
I_zz = (81/12)*a4 - 2 * ( (1/12)*a4 + a4 )
I_zz = (81/12)*a4 - 2 * (13/12)*a4
I_zz = (81 - 26) / 12 * a4
I_zz = 55 / 12 * a4 # I_zz = 55

# For I_ss, the distance 'd' is in the s-direction. d1 = a/2, d2 = -a/2
d_s = a/2
I_ss = I_main - 2 * (I_cutout_c + Area_cutout * d_s**2)
I_ss = (81/12)*a4 - 2 * ( (1/12)*a4 + a*a*(a/2)**2 )
I_ss = (81/12)*a4 - 2 * ( (1/12)*a4 + (1/4)*a4 )
I_ss = (81/12)*a4 - 2 * ( (1/12 + 3/12)*a4 )
I_ss = (81/12)*a4 - 2 * (4/12)*a4
I_ss = (81 - 8) / 12 * a4
I_ss = 73 / 12 * a4 # I_ss = 73

# Step 3: Calculate L and q_0 using the given data and calculated I_zz, I_ss
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

# Step 4: Calculate the force F using the derived formula F = (13 * q_0 * L) / 270
numerator_factor = 13
denominator = 270
F = (numerator_factor * q0 * L) / denominator

# Step 5: Print the results including the numbers in the final equation
print("The relationship between F, q_0, and L is: F = (13 * q_0 * L) / 270")
print(f"\nFirst, we calculate the parameters from the given geometry:")
print(f"I_zz = {I_zz}")
print(f"I_ss = {I_ss}")
print(f"L = (30 * {I_zz}) / 73 = {L}")
print(f"q_0 = (9 * {I_ss}) / 55 = {q0}")
print("\nSubstituting these values into the force equation:")
print(f"F = ({numerator_factor} * {q0} * {L}) / {denominator}")
print(f"F = {F}")
