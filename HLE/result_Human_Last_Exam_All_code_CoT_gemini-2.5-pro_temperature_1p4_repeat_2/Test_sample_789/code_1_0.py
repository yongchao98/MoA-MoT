import math

# --- Calculation Derivation ---
# Define the initial constants as 5-bit fractions
rho_n, rho_d = 9, 10
r_n, r_d = 1, 2
four_thirds_n, four_thirds_d = 4, 3
pi_n, pi_d = 22, 7

# Calculate r^3
r3_n, r3_d = r_n**3, r_d**3

print("--- Calculation Derivation ---")
print("The formula for the mass of a sphere is: m = rho * (4/3) * pi * r^3")
print("Using 5-bit fractional approximations, the initial equation is:")
print(f"m = ({rho_n}/{rho_d}) * ({four_thirds_n}/{four_thirds_d}) * ({pi_n}/{pi_d}) * ({r3_n}/{r3_d})")

print("\nTo keep all intermediate numbers within the 5-bit limit (0-31), we calculate in steps:")

# Step 1: Combine (rho) and (4/3)
# (9/10) * (4/3) = 36/30. This is invalid.
# We must simplify by dividing the numerator and denominator by their greatest common divisor (6).
term1_final_n = 6
term1_final_d = 5
print(f"Step 1: Combine rho and 4/3 -> ({rho_n}/{rho_d}) * ({four_thirds_n}/{four_thirds_d}). The unsimplified result 36/30 is invalid. The simplified result is {term1_final_n}/{term1_final_d}.")

# Step 2: Combine (pi) and (r^3)
# (22/7) * (1/8) = 22/56. This is invalid.
# We must simplify by dividing the numerator and denominator by their greatest common divisor (2).
term2_final_n = 11
term2_final_d = 28
print(f"Step 2: Combine pi and r^3 -> ({pi_n}/{pi_d}) * ({r3_n}/{r3_d}). The unsimplified result 22/56 is invalid. The simplified result is {term2_final_n}/{term2_final_d}.")

# Step 3: Multiply the intermediate results
# (6/5) * (11/28) -> 66/140. Both numbers are > 31.
# The simplified fraction is 33/70, which is also invalid as 33 and 70 are > 31.
print(f"\nStep 3: Multiply the intermediate results -> ({term1_final_n}/{term1_final_d}) * ({term2_final_n}/{term2_final_d}).")
print("This would result in 66/140, which cannot be represented. According to the rules, we must replace one of the fractions with a less precise one.")

# We find the best approximation for 11/28 (approx 0.3928) that allows the calculation to complete.
# The fraction 7/18 (approx 0.3888) is a close and suitable choice.
approx_n, approx_d = 7, 18
print(f"We replace ({term2_final_n}/{term2_final_d}) with a close approximation: ({approx_n}/{approx_d}).")

# New Step 3: (6/5) * (7/18)
# We simplify this expression before multiplying to stay within limits:
# (6 * 7) / (5 * 18) = (6 * 7) / (5 * 3 * 6)
# Cancelling the 6 gives 7 / (5 * 3) = 7/15.
final_n, final_d = 7, 15
print(f"The new calculation is ({term1_final_n}/{term1_final_d}) * ({approx_n}/{approx_d}), which simplifies to {final_n}/{final_d}.")
print(f"This is a valid 5-bit fraction. Our calculated mass is {final_n}/{final_d} kg.")

print("\n--- Final Equation ---")
print("The final calculation performed is:")
print(f"({term1_final_n} / {term1_final_d}) * ({approx_n} / {approx_d}) = {final_n} / {final_d}")

# --- Error Analysis ---
# Calculate the true value using high precision for pi
true_mass = 0.9 * (4.0/3.0) * math.pi * (0.5**3)
calculated_mass = final_n / final_d
abs_error = abs(true_mass - calculated_mass)
rounded_error = round(abs_error, 3)

print("\n--- Error Analysis ---")
print(f"True mass (high precision): {true_mass}")
print(f"Calculated mass (Titan 5-bit): {calculated_mass}")
print(f"Absolute error: |{true_mass} - {calculated_mass}| = {abs_error}")
print(f"The smallest absolute error, rounded to 0.001, is: {rounded_error}")

print(f"\n<<<{rounded_error}>>>")