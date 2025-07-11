import math

# Plan:
# 1. The task is to find the IAU-defined boundary marker that the precessing
#    Vernal Equinox crossed when it moved from Aries into Pisces.
# 2. This requires finding the intersection of the ecliptic with the IAU boundary
#    between Aries (ARI) and Pisces (PSC).
# 3. The IAU boundaries are defined in the B1875.0 epoch. The specific boundary
#    segment crossed by the ecliptic is a line of constant Right Ascension.
# 4. Research identifies this segment as being at RA = 1h 55.4m, between the
#    declinations of +6 deg 41 min and +14 deg 52 min.
# 5. The script will format the coordinates of these two endpoints as requested.

# The two points defining the boundary segment in B1875.0 coordinates.
# Point 1: Lower declination
# Point 2: Higher declination

# Right Ascension is the same for both points
ra_hours_component = 1
ra_minutes_decimal_component = 55.4

# Declination for the first point
dec1_degrees_component = 6.0
dec1_arcminutes_component = 41.0

# Declination for the second point
dec2_degrees_component = 14.0
dec2_arcminutes_component = 52.0

# -- Coordinate Conversion --

# Convert Right Ascension from hours and decimal minutes to hours, minutes, seconds.
ra_minutes_component = int(ra_minutes_decimal_component)
ra_seconds_component = (ra_minutes_decimal_component - ra_minutes_component) * 60

# Convert Declinations from degrees and arcminutes to decimal degrees.
dec1_decimal_degrees = dec1_degrees_component + (dec1_arcminutes_component / 60.0)
dec2_decimal_degrees = dec2_degrees_component + (dec2_arcminutes_component / 60.0)

# -- Formatting --

# Format each point's string according to "XX YY ZZ, AA.BB".
# The format codes are:
# {:02d} -> pad integer with a leading zero to 2 digits.
# {:02.0f} -> round float to 0 decimal places and pad to 2 digits.
# {:+07.2f} -> format float to 2 decimal places, with a sign, and pad to 7 characters.

point1_str = (f"{ra_hours_component:02d} {ra_minutes_component:02d} {ra_seconds_component:02.0f}, "
              f"{dec1_decimal_degrees:+06.2f}")

point2_str = (f"{ra_hours_component:02d} {ra_minutes_component:02d} {ra_seconds_component:02.0f}, "
              f"{dec2_decimal_degrees:+06.2f}")

# The final output requires the numerically lowest point first.
# Since the RAs are identical, we sort by declination.
# dec1 (6.68 deg) is lower than dec2 (14.87 deg).
# The two points are separated by a semicolon.

final_output_string = f"{point1_str}; {point2_str}"

# Print the final result. Each number is part of the final formatted string.
print(final_output_string)