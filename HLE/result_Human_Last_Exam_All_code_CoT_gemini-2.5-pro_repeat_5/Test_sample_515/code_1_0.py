import math

# Step 1: Define the constants and parameters of the problem.
# The starting point is x0 = (0,1). The target point is z = (3600,0).
# The probability p is approximated by a(x0) / a(z).

# The value of the potential kernel at a nearest neighbor of the origin is known.
a_x0 = 1.0
print(f"The potential kernel at the starting point (0,1) is a((0,1)) = {a_x0}")

# The distance of the target point z from the origin.
z_norm = 3600.0

# Constants for the asymptotic formula of the potential kernel a(z).
# Euler-Mascheroni constant
gamma = 0.5772156649
# Natural logarithm of 8
ln8 = math.log(8)

# Step 2: Calculate the potential kernel at the target point z using the asymptotic formula.
# a(z) ~ (2/pi) * (ln(|z|) + gamma + ln(8)) for large |z|.
ln_z = math.log(z_norm)
a_z = (2 / math.pi) * (ln_z + gamma + ln8)

print(f"The distance of the target point from the origin is |z| = {z_norm}")
print(f"The natural logarithm of the distance is ln(|z|) = {ln_z}")
print(f"The Euler-Mascheroni constant is gamma = {gamma}")
print(f"The natural logarithm of 8 is ln(8) = {ln8}")
print(f"The potential kernel at the target point (3600,0) is approximately a(z) = (2/pi) * ({ln_z} + {gamma} + {ln8}) = {a_z}")

# Step 3: Calculate the probability.
p = a_x0 / a_z
print(f"\nThe probability is p = a(x0) / a(z) = {a_x0} / {a_z} = {p}")

# Step 4: Give the approximate answer with two significant digits.
p_approx = float(f"{p:.2g}")
print(f"\nThe approximate answer with two significant digits is {p_approx}")
