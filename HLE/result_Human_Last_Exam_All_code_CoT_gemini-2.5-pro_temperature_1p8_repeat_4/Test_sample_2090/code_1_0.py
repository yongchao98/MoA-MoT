# The apparent position of the Sun at the Vernal Equinox is at RA=0, Dec=0.
# Due to precession, this point moves. It crossed from Aries into Pisces.
# The IAU boundary between these constellations at the celestial equator is a line
# of constant Right Ascension defined for the B1875.0 epoch.

# This boundary is located at Right Ascension 0 hours.
# It is a line segment defined by its endpoints, which are junctions with
# other constellations.

# Northern endpoint (at the boundary with Andromeda)
ra_h_north, ra_m_north, ra_s_north = 0, 0, 0
dec_north = 6.0

# Southern endpoint (at the boundary with Cetus)
ra_h_south, ra_m_south, ra_s_south = 0, 0, 0
dec_south = -9.0

# The request asks for the "numerically lowest point" to be listed first.
# -9.0 is numerically lower than 6.0, so we use the southern point first.
point1_ra_str = f"{ra_h_south:02d} {ra_m_south:02d} {ra_s_south:02d}"
point1_dec_str = f"{dec_south:.2f}"
point1_full_str = f"{point1_ra_str}, {point1_dec_str}"

point2_ra_str = f"{ra_h_north:02d} {ra_m_north:02d} {ra_s_north:02d}"
point2_dec_str = f"{dec_north:.2f}"
point2_full_str = f"{point2_ra_str}, {point2_dec_str}"

# Output each number/component of the points for clarity, as requested.
print(f"The relevant boundary is a line segment at Right Ascension {ra_h_south}h {ra_m_south}m {ra_s_south}s.")
print(f"The southern endpoint has a Declination of {dec_south:.2f} degrees.")
print(f"The northern endpoint has a Declination of {dec_north:.2f} degrees.")
print("\nFormatted boundary marker points:")

# Combine the two points into the final answer string, separated by a semicolon.
final_answer_string = f"{point1_full_str}; {point2_full_str}"

print(final_answer_string)
