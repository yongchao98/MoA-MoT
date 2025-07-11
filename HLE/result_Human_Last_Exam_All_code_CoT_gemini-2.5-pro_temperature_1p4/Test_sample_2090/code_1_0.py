# The coordinates for the boundary between Aries and Pisces
# are based on the IAU B1875.0 standard definition.
# The boundary is a line of constant Right Ascension at 0h.
# The specific segment crossing the celestial equator runs from Dec -15 to +5 degrees.

# Define the coordinates of the two endpoints of the boundary line segment.
# Point 1 (the "lower" point based on declination)
ra_h_1, ra_m_1, ra_s_1 = 0, 0, 0
dec_1 = -15.0

# Point 2 (the "higher" point)
ra_h_2, ra_m_2, ra_s_2 = 0, 0, 0
dec_2 = 5.0

# Format the points according to the "XX YY ZZ, AA.BB" specification.
# The zfill(2) function pads the numbers with a leading zero if they are less than 10.
point_1_str = f"{str(ra_h_1).zfill(2)} {str(ra_m_1).zfill(2)} {str(ra_s_1).zfill(2)}, {dec_1:.2f}"
point_2_str = f"{str(ra_h_2).zfill(2)} {str(ra_m_2).zfill(2)} {str(ra_s_2).zfill(2)}, {dec_2:.2f}"

# Combine the two points with a semicolon, with the numerically lowest point first.
# Since RAs are identical, we order by declination. -15.00 comes before 5.00.
final_answer = f"{point_1_str}; {point_2_str}"

# Print the final result.
print(final_answer)
