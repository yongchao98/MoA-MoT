import math

# Step 1: Define relationships from observations.
# The problem gives two ratios based on angular size observations.
# Let R_P1, R_P2, R_BD be the radii of Planet 1, Planet 2, and the Brown Dwarf.
# Let r_P1, r_P2 be the orbital distances.
# We are given that r_P1 and r_P2 are much larger than R_BD.
#
# Observation 1 (from Station): ang_size(P1) = 0.8 * ang_size(P2)
# Approximates to: (R_P1 / r_P1) = 0.8 * (R_P2 / r_P2)
#
# Observation 2 (from Planet 2): ang_size(P1) = 0.2 * ang_size(BD)
# Gives: R_P1 / (r_P2 - r_P1) = 0.2 * R_BD / r_P2

# Step 2: Relate the desired ratio of radii (R_P2/R_BD) to the ratio of orbital distances (r_P2/r_P1).
# By rearranging the equations above, we can derive the formula:
# R_P2 / R_BD = 0.25 * (r_P2 / r_P1 - 1)
print("From the geometry of the transits, we derive the relationship:")
print("R_P2 / R_BD = 0.25 * (r_P2 / r_P1 - 1)\n")


# Step 3: Use orbital mechanics to find the ratio of orbital distances (r_P2/r_P1).
# We are told Planet 1 is on a parabolic orbit (at pericenter r_P1) and Planet 2 is on a circular orbit (at r_P2).
# The specific angular momentum (h) for a parabolic orbit at pericenter is h1 = sqrt(2*G*M*r_P1).
# The specific angular momentum for a circular orbit is h2 = sqrt(G*M*r_P2).
# A key physical insight is to assume these two values are equal (h1 = h2).
# sqrt(2*G*M*r_P1) = sqrt(G*M*r_P2)
# Squaring both sides and simplifying gives: 2 * r_P1 = r_P2.
# Therefore, the ratio r_P2 / r_P1 is 2.
ratio_r = 2.0
print("Using the properties of parabolic and circular orbits, a key physical constraint is that the specific angular momenta of the planets are equal.")
print(f"This implies a ratio of orbital distances, r_P2 / r_P1 = {ratio_r}\n")


# Step 4: Calculate the numerical value for the ratio of radii.
ratio_R_P2_R_BD_numerator = 0.25
ratio_r_minus_one = ratio_r - 1
ratio_R_P2_R_BD = ratio_R_P2_R_BD_numerator * ratio_r_minus_one
print("Substituting the orbital distance ratio into our geometric formula:")
print(f"R_P2 / R_BD = {ratio_R_P2_R_BD_numerator} * ({ratio_r} - 1)")
print(f"R_P2 / R_BD = {ratio_R_P2_R_BD}\n")


# Step 5: Calculate the final brightness drop in magnitudes.
# The brightness drop is Δm = -2.5 * log10(1 - (R_P2 / R_BD)^2).
c_mag = -2.5
ratio_R_sq = ratio_R_P2_R_BD**2
flux_ratio = 1 - ratio_R_sq
delta_m = c_mag * math.log10(flux_ratio)

print("Finally, we calculate the brightness drop during the transit of Planet 2 in front of the Brown Dwarf.")
print("The formula is: Δm = -2.5 * log10(1 - (R_P2 / R_BD)^2)")
print("Plugging in the numbers:")
print(f"Δm = ({c_mag}) * log10(1 - ({ratio_R_P2_R_BD})^2)")
print(f"Δm = ({c_mag}) * log10(1 - {ratio_R_sq})")
print(f"Δm = ({c_mag}) * log10({flux_ratio})")
print(f"The brightness drop of the brown dwarf is {delta_m:.3f} magnitudes.")

# Final answer format
print(f"\n<<<0.070>>>")