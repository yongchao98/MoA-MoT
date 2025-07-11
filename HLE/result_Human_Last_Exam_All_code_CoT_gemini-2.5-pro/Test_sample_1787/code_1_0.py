import math

# Step 1: Define the given ratios from the problem statement.
# ratio_p1_p2: Angular size of Planet 1 is 80% of Planet 2 from the station.
ratio_p1_p2 = 0.80
# ratio_p1_star: Angular size of Planet 1 is 20% of the brown dwarf from Planet 2.
ratio_p1_star = 0.20

# Step 2: Establish the relationship between the orbital radii.
# Based on orbital mechanics, we assume the specific angular momentum of the two
# planets are equal (h1 = h2).
# h1 = sqrt(2*G*M*r1) for a parabolic orbit at pericenter r1.
# h2 = sqrt(G*M*r2) for a circular orbit at radius r2.
# Equating them (h1=h2) gives sqrt(2*G*M*r1) = sqrt(G*M*r2), which simplifies to r2 = 2*r1.
r2_over_r1 = 2.0

# Step 3: Solve for the ratio of Planet 2's radius to the Brown Dwarf's radius (R2/Rs).
# From observation 1: R1/r1 = ratio_p1_p2 * (R2/r2)
# From observation 2: R1/(r2-r1) = ratio_p1_star * (Rs/r2)
# We can express R2/Rs in terms of the known ratios:
# R2/Rs = (ratio_p1_star / ratio_p1_p2) * (r2 / (r2 - r1))
# R2/Rs = (ratio_p1_star / ratio_p1_p2) * (r2/r1 / (r2/r1 - 1))
ratio_R2_Rs = (ratio_p1_star / ratio_p1_p2) * (r2_over_r1 / (r2_over_r1 - 1.0))

# Step 4: Calculate the squared ratio, which corresponds to the ratio of the disks' areas.
ratio_area = ratio_R2_Rs**2

# Step 5: Calculate the flux ratio during the transit.
flux_ratio = 1 - ratio_area

# Step 6: Calculate the brightness drop in bolometric magnitudes.
delta_magnitude = -2.5 * math.log10(flux_ratio)

# Step 7: Print the final equation with the calculated values and the result.
print("The brightness drop (Δm) is calculated using the formula:")
print(f"Δm = -2.5 * log10(1 - (R_planet2 / R_star)^2)")
print("\nFirst, we find the ratio of the radii, R_planet2 / R_star:")
print(f"R_planet2 / R_star = {ratio_R2_Rs:.4f}")
print("\nPlugging the values into the formula:")
print(f"Δm = -2.5 * log10(1 - ({ratio_R2_Rs:.4f})^2)")
print(f"Δm = -2.5 * log10(1 - {ratio_area:.4f})")
print(f"Δm = -2.5 * log10({flux_ratio:.4f})")
print(f"\nThe calculated brightness drop is: {delta_magnitude:.3f} magnitudes.")

# Final answer in the specified format
print(f"\n<<<{delta_magnitude:.3f}>>>")