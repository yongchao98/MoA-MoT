import math

# Step 1: Define the constants and the starting/target points.
x0 = (0, 1)
z = (3600, 0)

# Step 2: Calculate the potential kernel a(x) for the start and target points.

# For a neighbor of the origin, a(x0) is 1 by standard normalization.
a_x0 = 1.0
print(f"The potential kernel at the starting point (0,1) is a(0,1) = {a_x0}")

# For the target point z, |z| = 3600, which is large.
# We use the asymptotic formula: a(z) ≈ (2/π) * ln|z| + C
# where C = (2γ + ln(8))/π.
gamma = 0.5772156649  # Euler-Mascheroni constant
pi = math.pi
ln8 = math.log(8)
z_norm = 3600.0

# Calculate the constant C
C = (2 * gamma + ln8) / pi
# Calculate a(z)
ln_z_norm = math.log(z_norm)
a_z = (2 / pi) * ln_z_norm + C

print(f"The potential kernel at the target point (3600,0) is a(3600,0) ≈ (2/π) * ln(3600) + C")
print(f"where C = (2γ + ln(8))/π ≈ {C:.4f}")
print(f"a(3600,0) ≈ (2/{pi:.4f}) * {ln_z_norm:.4f} + {C:.4f} ≈ {a_z:.4f}")

# Step 3: Calculate the probability.
# p ≈ a(x0) / a(z)
prob = a_x0 / a_z

print(f"\nThe probability is given by the formula p ≈ a(0,1) / a(3600,0)")
print(f"p ≈ {a_x0} / {a_z:.4f}")
print(f"p ≈ {prob:.4f}")

# Final answer with two significant digits
prob_approx = float(f"{prob:.2g}")
print(f"\nThe approximate answer with two significant digits is {prob_approx}")
