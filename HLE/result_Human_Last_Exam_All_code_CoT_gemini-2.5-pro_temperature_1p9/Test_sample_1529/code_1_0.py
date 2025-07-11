# The given curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1 is a hyperelliptic curve.
# It is known to be birationally equivalent to the elliptic curve Y^2 = X^3 - 4X + 16.
# This is its minimal Weierstrass model. We will calculate the discriminant of this model.

# Coefficients of the minimal Weierstrass equation y^2 = x^3 + Ax + B
A = -4
B = 16

# The formula for the discriminant is Delta = -16 * (4*A^3 + 27*B^2)

# Perform the calculation
val_A_cubed = 4 * (A**3)
val_B_squared = 27 * (B**2)
discriminant = -16 * (val_A_cubed + val_B_squared)

# Print the information for the user
print("The original curve is birationally equivalent to the elliptic curve with the equation:")
print(f"y^2 = x^3 + ({A})x + ({B})")
print("\nThe minimal discriminant is calculated using the formula Delta = -16 * (4*A^3 + 27*B^2).")
print("\nCalculation steps:")
print(f"A = {A}")
print(f"B = {B}")
print(f"4 * A^3 = 4 * ({A})^3 = 4 * {A**3} = {val_A_cubed}")
print(f"27 * B^2 = 27 * ({B})^2 = 27 * {B**2} = {val_B_squared}")
print(f"Delta = -16 * ({val_A_cubed} + {val_B_squared})")
print(f"Delta = -16 * ({val_A_cubed + val_B_squared})")
print(f"Delta = {discriminant}")

print("\nThe final minimal discriminant is:")
print(discriminant)