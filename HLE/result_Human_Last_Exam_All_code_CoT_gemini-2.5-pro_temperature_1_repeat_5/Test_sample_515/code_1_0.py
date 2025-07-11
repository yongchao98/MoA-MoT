import math

# Step 6 & 7: Calculate the values and the final probability.

# The starting point is x0 = (0,1).
# The target set A is the set of neighbors of z = (3600,0).
# The probability p is approximated by a(x0) / a(z).

# For the potential kernel a(x) normalized such that a(0)=0, we have a(1,0) = 1.
# By symmetry, a(0,1) is the same as a(1,0).
a_x0 = 1.0
print(f"The potential kernel at the starting point (0,1) is a(0,1) = {a_x0}")

# For the distant point z=(3600,0), we use the asymptotic formula a(z) ~ (2/pi) * log(|z|).
# The distance |z| is 3600.
z_norm = 3600.0
a_z = (2 / math.pi) * math.log(z_norm)
print(f"The potential kernel at the target point (3600,0) is approximately a(3600,0) = (2/pi) * log(3600) = {a_z:.4f}")

# The probability is the ratio of these two values.
probability = a_x0 / a_z

print(f"\nThe probability is approximately given by the equation: p = a(0,1) / a(3600,0)")
print(f"p = {a_x0} / {a_z:.4f}")
print(f"p ≈ {probability:.4f}")

# Give the approximate answer with two significant digits.
print(f"\nThe approximate answer with two significant digits is: {probability:.2g}")