import math

# --- Problem Setup ---
# In frame F, for any pair i, j, cos(theta_ij) = -1/3.
one_minus_cos_theta_F = 1 - (-1/3)

# In frame F', we are given the following angles and their cosines:
# theta'_12 = pi/2 => cos = 0
one_minus_cos_theta_prime_12 = 1 - 0
# theta'_13 = 3pi/4 => cos = -sqrt(2)/2
one_minus_cos_theta_prime_13 = 1 - (-math.sqrt(2)/2)
# theta'_23 = 3pi/4 => cos = -sqrt(2)/2
one_minus_cos_theta_prime_23 = 1 - (-math.sqrt(2)/2)

print("This script calculates the value of (1 - cos(theta'_14)) / (1 - cos(theta'_34)).\n")
print("This ratio is equivalent to (E'_3/E3) / (E'_1/E1), which we can find using the given angles.\n")

# --- Solving for Energy Ratios ---
# We have a system of equations for the energy ratios r_i = E'_i/E_i:
# r1 * r2 = one_minus_cos_theta_F / one_minus_cos_theta_prime_12
# r1 * r3 = one_minus_cos_theta_F / one_minus_cos_theta_prime_13
# r2 * r3 = one_minus_cos_theta_F / one_minus_cos_theta_prime_23

# Right-hand side of the equations
RHS1 = one_minus_cos_theta_F / one_minus_cos_theta_prime_12
RHS2 = one_minus_cos_theta_F / one_minus_cos_theta_prime_13
RHS3 = one_minus_cos_theta_F / one_minus_cos_theta_prime_23

# From RHS2 = RHS3, we deduce r1 = r2.
# Substitute r1=r2 into the first equation: r1^2 = RHS1.
r1_squared = RHS1
r1 = math.sqrt(r1_squared)

# Now find r3 using the second equation: r1 * r3 = RHS2.
r3 = RHS2 / r1

print("--- Calculation Steps ---")
print(f"1. From theta'_12 = pi/2, we find the ratio E'_1/E_1:")
print(f"   (E'_1/E_1)^2 = (1 - cos(theta_12)) / (1 - cos(theta'_12)) = {one_minus_cos_theta_F:.4f} / {one_minus_cos_theta_prime_12:.4f} = {r1_squared:.4f}")
print(f"   So, E'_1/E_1 = sqrt({r1_squared:.4f}) = {r1:.4f}\n")

print(f"2. From theta'_13 = 3pi/4, we find the ratio E'_3/E_3:")
print(f"   (E'_1/E_1) * (E'_3/E_3) = (1 - cos(theta_13)) / (1 - cos(theta'_13)) = {one_minus_cos_theta_F:.4f} / {one_minus_cos_theta_prime_13:.4f} = {RHS2:.4f}")
print(f"   So, E'_3/E_3 = {RHS2:.4f} / {r1:.4f} = {r3:.4f}\n")


# --- Final Result ---
final_value = r3 / r1

print("--- Final Result ---")
print("The value of the expression is given by the equation:")
print("(1 - cos(\u03B8'_{14})) / (1 - cos(\u03B8'_{34})) = (E'_3/E_3) / (E'_1/E_1)")
print(f"                                   = {r3:.5f} / {r1:.5f}")
print(f"                                   = {final_value:.5f}\n")

# For verification, the exact answer is 2 - sqrt(2)
exact_answer = 2 - math.sqrt(2)
print(f"The exact value is 2 - sqrt(2), which is approximately {exact_answer:.5f}.")
