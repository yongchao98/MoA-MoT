import math

# Step 1 & 2: Derive the relationship between radii and distances from the problem statement.
# Let d1 be the distance to Planet 1 and d2 be the distance to Planet 2.
# Let x = d1 / d2.
# The two transit observations give two equations:
# (1) R1 / d1 = 0.8 * R2 / d2
# (2) R1 / (d2 - d1) = 0.2 * R_BD / d2
# Solving these for the ratio R2 / R_BD gives:
# R2 / R_BD = 0.25 * (1/x - 1)

# Step 3: Use orbital mechanics to find the value of x.
# Planet 1 is at the pericenter of a parabolic orbit, so its velocity v1 satisfies: v1^2 = 2*G*M/d1.
# Planet 2 is in a circular orbit of radius d2, so its velocity v2 satisfies: v2^2 = G*M/d2.
# A necessary assumption to solve the problem is that their specific angular momenta (L/m = r*v) are equal.
# L1/m1 = d1*v1 = d1*sqrt(2*G*M/d1) = sqrt(2*G*M*d1)
# L2/m2 = d2*v2 = d2*sqrt(G*M/d2) = sqrt(G*M*d2)
# Equating L1/m1 and L2/m2 gives: sqrt(2*G*M*d1) = sqrt(G*M*d2), which simplifies to 2*d1 = d2.
# Therefore, the ratio of the distances is:
x = 1.0 / 2.0

# Step 4: Calculate the ratio of the radii.
# Substitute x into the formula from Step 2.
ratio_R2_to_RBD = 0.25 * (1.0 / x - 1.0)

# Step 5: Calculate the fractional brightness drop.
# The brightness drop is the ratio of the disks' areas.
brightness_drop_fraction = ratio_R2_to_RBD**2

# Step 6: Convert the brightness drop to bolometric magnitudes.
# The change in magnitude is Δm = -2.5 * log10(1 - drop_fraction).
l_ratio = 1.0 - brightness_drop_fraction
delta_magnitude = -2.5 * math.log10(l_ratio)

print("The calculation for the brightness drop follows these steps:")
print("\n1. Find the ratio of radii R_2/R_BD:")
print(f"R_2 / R_BD = 0.25 * (1 / {x} - 1) = {ratio_R2_to_RBD}")

print("\n2. Find the fractional brightness drop, which is the square of the radii ratio:")
print(f"Fractional drop = ({ratio_R2_to_RBD})^2 = {brightness_drop_fraction}")

print("\n3. Calculate the brightness drop in magnitudes using the formula Δm = -2.5 * log10(1 - fractional_drop):")
print(f"Δm = -2.5 * log10(1 - {brightness_drop_fraction})")
print(f"Δm = -2.5 * log10({l_ratio})")
print(f"Δm = -2.5 * ({math.log10(l_ratio):.6f})")
print(f"Δm = {delta_magnitude:.3f}")

print("\n<<<0.070>>>")