import numpy as np

# The problem simplifies to calculating the value of the expression:
# (8/15) * pi^8 + (1/3) * pi^2 - (1/2) * pi + 1
# This is derived by simplifying the denominator and splitting the integral into four parts,
# each corresponding to a known mathematical constant or formula.

# Define the coefficients of the final expression
final_coeff_pi8_num = 8
final_coeff_pi8_den = 15

final_coeff_pi2_num = 1
final_coeff_pi2_den = 3

final_coeff_pi_num = -1
final_coeff_pi_den = 2

final_const = 1

# Calculate the final numerical value
final_value = (final_coeff_pi8_num / final_coeff_pi8_den) * np.pi**8 \
            + (final_coeff_pi2_num / final_coeff_pi2_den) * np.pi**2 \
            + (final_coeff_pi_num / final_coeff_pi_den) * np.pi \
            + final_const

print("The integral evaluates to the following expression:")
print(f"I = ({final_coeff_pi8_num}/{final_coeff_pi8_den})*pi^8 + ({final_coeff_pi2_num}/{final_coeff_pi2_den})*pi^2 + ({final_coeff_pi_num}/{final_coeff_pi_den})*pi + {final_const}")
print("\nExplicitly showing each number in the final equation:")
print(f"Value = ({final_coeff_pi8_num} / {final_coeff_pi8_den}) * pi^8 + ({final_coeff_pi2_num} / {final_coeff_pi2_den}) * pi^2 - ({abs(final_coeff_pi_num)} / {final_coeff_pi_den}) * pi + {final_const}")

print(f"\nThe numerical value of the integral is approximately: {final_value}")