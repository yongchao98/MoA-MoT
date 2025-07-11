import math

# Angular size ratio from the station: Planet 1 vs Planet 2
k_S = 0.8
# Angular size ratio from Planet 2: Planet 1 vs Brown Dwarf
k_P2 = 0.2

# The relationship between the ratio of radii (Q = R_P2/R_BD) and the ratio of orbital distances (x = r_P1/r_P2)
# is derived as Q = (k_P2 / k_S) * (1-x)/x
ratio_k = k_P2 / k_S
print(f"The ratio of the observational constants k_P2/k_S is {k_P2:.1f}/{k_S:.1f} = {ratio_k:.2f}")

# From the orbital mechanics (assuming equal specific angular momentum for the parabolic and circular orbits),
# we found the ratio of the orbital distances: r_P2 = 2 * r_P1.
# Therefore, x = r_P1 / r_P2 = 1/2.
x = 0.5
print(f"The ratio of the orbital distances, x = r_P1/r_P2, is found to be {x}")

# Now we can calculate the ratio of the radii of Planet 2 and the Brown Dwarf.
ratio_radii = ratio_k * (1 - x) / x
print(f"The ratio of the radii, R_P2/R_BD = {ratio_k:.2f} * (1 - {x}) / {x} = {ratio_radii:.2f}")

# The fractional brightness drop (f) is the square of the radii ratio.
f = ratio_radii**2
print(f"The fractional brightness drop f = (R_P2/R_BD)^2 = {ratio_radii:.2f}^2 = {f:.4f}")

# The brightness drop in bolometric magnitudes (delta_m) is given by the formula:
# delta_m = -2.5 * log10(1 - f)
delta_m = -2.5 * math.log10(1 - f)

print("\nThe final equation for the brightness drop in magnitudes is:")
print(f"Δm = -2.5 * log10(1 - {ratio_radii:.2f}^2)")
print(f"Δm = -2.5 * log10(1 - {f:.4f})")
print(f"Δm = -2.5 * log10({1-f:.4f})")
print(f"Δm = {delta_m:.3f}")

# Final Answer
print("\nThe brightness drop of the brown dwarf is:")
print(f"{delta_m:.3f} magnitudes")
<<<0.070>>>