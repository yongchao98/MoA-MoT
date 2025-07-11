# Plan:
# 1. Define the coordinates of the two vertices based on the IAU 1875.0 boundary definitions.
# 2. Convert the Right Ascension from hours and decimal minutes into hours, minutes, and seconds.
# 3. Convert the Declinations from degrees and arcminutes into decimal degrees.
# 4. Format the two points into the required "XX YY ZZ, AA.BB" string format.
# 5. Arrange the points with the numerically lowest declination first.
# 6. Print the final combined string.

# Step 1: Define the coordinates of the two vertices.
# These values are from the IAU's 1930 definition of constellation boundaries.
# The relevant boundary is a line of constant Right Ascension.
ra_h = 23
ra_m_decimal = 58.4

# The first point's declination (southern vertex, Pisces/Cetus border)
dec1_deg = -6
dec1_arcmin = 0

# The second point's declination (northern vertex, Pisces/Aries border)
dec2_deg = 5
dec2_arcmin = 30

# Step 2: Convert Right Ascension to H:M:S format.
ra_m_integer = int(ra_m_decimal)
ra_s_float = (ra_m_decimal - ra_m_integer) * 60
ra_s_int = int(round(ra_s_float))

# Step 3: Convert Declinations to decimal degree format.
dec1_decimal = float(dec1_deg) + float(dec1_arcmin) / 60.0
dec2_decimal = float(dec2_deg) + float(dec2_arcmin) / 60.0

# Step 4: Format the two points into strings.
# The format is specified as "XX YY ZZ, AA.BB".
# The point with the lower declination (-6.00) is defined first.
point1_str = f"{ra_h:02d} {ra_m_integer:02d} {ra_s_int:02d}, {dec1_decimal:.2f}"

# The point with the higher declination (5.50) is defined second.
point2_str = f"{ra_h:02d} {ra_m_integer:02d} {ra_s_int:02d}, {dec2_decimal:.2f}"

# Step 5 & 6: Combine the two point strings, separated by a semicolon, and print.
# The problem asks for the numerically lowest point first, which is the one with declination -6.00.
final_output = f"{point1_str}; {point2_str}"
print(final_output)
