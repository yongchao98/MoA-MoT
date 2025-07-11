import numpy as np

# Step 1: Coefficients for Phi(z) = a*z + b
# From 8a + b = 5 and 2a = 1
a = 1/2
b = 1

# Step 2: The sum of integrals I_1 + I_2 simplifies to 11/3
I_sum = 11/3

# Step 3: Calculate the total energy E_total = 1 + 0.5 * (I_1 + I_2)
# We use the simplified sum I_sum.
constant_term = 1
integral_contribution = 0.5 * I_sum

E_total = constant_term + integral_contribution

# Step 4: Output the equation with the calculated numbers
# E_total = 1 + (1/2) * (11/3) = 17/6
val_1 = 1
val_2 = 0.5
val_3 = I_sum
final_val = E_total

print(f"The minimum total heat energy is calculated as:")
print(f"E_total = {val_1} + {val_2} * ({val_3:.4f})")
print(f"E_total = {val_1} + {val_2 * val_3:.4f}")
print(f"E_total = {final_val:.4f}")

# Final Answer as a fraction
# E_total = 1 + 11/6 = 6/6 + 11/6 = 17/6
final_numerator = 17
final_denominator = 6
print(f"In fractional form, E_total = {final_numerator}/{final_denominator}")
<<<17/6>>>