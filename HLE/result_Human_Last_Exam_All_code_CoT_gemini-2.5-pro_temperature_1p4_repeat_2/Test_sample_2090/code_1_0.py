import math

def format_ra(hours):
    """Formats Right Ascension from decimal hours to HH MM SS."""
    h = int(hours)
    minutes_decimal = (hours - h) * 60
    m = int(minutes_decimal)
    s = int((minutes_decimal - m) * 60)
    return f"{h:02d} {m:02d} {s:02d}"

def format_dec(degrees):
    """Formats Declination to AA.BB."""
    return f"{degrees:.2f}"

# According to IAU constellation boundary data (epoch B1875.0), the
# boundary line separating Pisces and Aries is at Right Ascension 0h.
# The specific line segment that crosses the celestial equator (and thus the ecliptic at the equinox)
# is defined by two vertices.
# The northern vertex is at (RA=0h, Dec=+15.0 degrees).
# The southern vertex is the triple point between Aries, Pisces, and Cetus,
# located at (RA=0h, Dec=-0.9 degrees).

# Define the two points of the line segment
# Point 1: The southern vertex
ra1_h = 0.0
dec1_deg = -0.9

# Point 2: The northern vertex
ra2_h = 0.0
dec2_deg = 15.0

# The problem asks to list the numerically lowest point first.
# Since the RA is the same for both, we order by Declination.
point1_str = f"{format_ra(ra1_h)}, {format_dec(dec1_deg)}"
point2_str = f"{format_ra(ra2_h)}, {format_dec(dec2_deg)}"

# Combine and print the final answer string
final_answer = f"{point1_str}; {point2_str}"
print("The two points defining the boundary line are:")
print(final_answer)