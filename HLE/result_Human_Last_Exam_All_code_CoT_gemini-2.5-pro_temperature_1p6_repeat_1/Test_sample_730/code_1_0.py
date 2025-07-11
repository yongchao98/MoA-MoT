import math

# Step 1: Define the given values as fractions suitable for the parrot.
# Density (rho) = 0.9 = 9/10 kg/cm^3
rho_n, rho_d = 9, 10

# Radius (r) = 0.5 = 1/2 cm
r_n, r_d = 1, 2

# Step 2: Define the constants from the volume formula.
# Volume of a sphere = (4/3) * pi * r^3
vol_const_n, vol_const_d = 4, 3

# Step 3: Find a suitable fractional approximation for pi.
# The parrot can only handle integers up to 10.
# Let's test pi_approx = 10/3.
pi_approx_n, pi_approx_d = 10, 3

# Step 4: Write down the full calculation for the parrot.
# Mass = rho * (4/3) * pi * r^3
# The term r^3 is (1/2)^3 = 1/8.
r3_n, r3_d = r_n**3, r_d**3

# The integers involved in the expression are:
# For rho: 9, 10
# For volume constant: 4, 3
# For pi approximation: 10, 3
# For radius cubed: 1, 8 (since (1/2)^3 = 1/8)
# The largest integer is 10, which is acceptable.

# Step 5: Calculate the estimated mass to check the error.
mass_est_n = rho_n * vol_const_n * pi_approx_n * r3_n
mass_est_d = rho_d * vol_const_d * pi_approx_d * r3_d

# Simplify the resulting fraction
common_divisor = math.gcd(mass_est_n, mass_est_d)
final_n = mass_est_n // common_divisor
final_d = mass_est_d // common_divisor
estimated_mass = final_n / final_d

# Step 6: Calculate the exact mass for comparison.
exact_mass = (rho_n / rho_d) * (vol_const_n / vol_const_d) * math.pi * (r_n / r_d)**3
error = abs(estimated_mass - exact_mass) / exact_mass * 100

# Step 7: Print the instruction for the parrot.
print("Yes, the parrot can estimate the mass.")
print("Here are the instructions for the calculation:\n")
print("The formula is Mass = Density * (4/3) * pi * (Radius)^3")
print(f"Use Density = {rho_n}/{rho_d}, approximate pi as {pi_approx_n}/{pi_approx_d}, and use Radius = {r_n}/{r_d}.")
print("\nThe final calculation to perform is:")
# The problem asks to output each number in the final equation.
# So we present the full, unsimplified expression and its simple result.
print(f"({rho_n}/{rho_d}) * ({vol_const_n}/{vol_const_d}) * ({pi_approx_n}/{pi_approx_d}) * ({r_n}/{r_d})^3 = {final_n}/{final_d}")

print(f"\nThis estimation gives a mass of {estimated_mass:.3f} kg.")
print(f"The actual mass is {exact_mass:.3f} kg.")
print(f"The error is {error:.2f}%, which is within the 10% limit.")