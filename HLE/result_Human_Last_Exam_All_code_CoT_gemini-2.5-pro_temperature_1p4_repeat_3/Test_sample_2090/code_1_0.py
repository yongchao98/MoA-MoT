# The coordinates of the two points defining the boundary segment between
# Pisces and Aries that crosses the celestial equator.
# The data is based on the IAU's B1875.0 epoch definitions.

# Point 1 (Southern endpoint, numerically lower declination)
ra_h_1 = 23
ra_m_1 = 48
ra_s_1 = 30
dec_1_deg = -2.70

# Point 2 (Northern endpoint)
ra_h_2 = 23
ra_m_2 = 48
ra_s_2 = 30
dec_2_deg = 4.30

# The problem requires printing the final formatted string.
# The format is "RA(h m s), Dec(deg.dd); RA(h m s), Dec(deg.dd)"
# The two points are separated by a semicolon.
# The point with the numerically lowest value appears first.
# Since the Right Ascension is the same, we order by Declination.

print(f"{ra_h_1} {ra_m_1} {ra_s_1}, {dec_1_deg:.2f}; {ra_h_2} {ra_m_2} {ra_s_2}, {dec_2_deg:.2f}")