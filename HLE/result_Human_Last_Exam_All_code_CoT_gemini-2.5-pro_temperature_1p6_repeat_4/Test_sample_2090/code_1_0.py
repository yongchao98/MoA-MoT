import math

def format_ra_hms(ra_hours):
    """Converts RA from decimal hours to HH MM SS format."""
    h = int(ra_hours)
    minutes_decimal = (ra_hours - h) * 60
    m = int(minutes_decimal)
    s = int((minutes_decimal - m) * 60)
    return f"{h:02d}", f"{m:02d}", f"{s:02d}"

def format_dec_deg(dec_degrees):
    """Formats declination to two decimal places."""
    return f"{dec_degrees:.2f}"

def solve_boundary_crossing():
    """
    Identifies the IAU B1875.0 boundary between Aries and Pisces crossed by the
    ecliptic and prints its coordinates.
    """
    # Step 1: Define the coordinates of the Aries/Pisces boundary segment
    # These are from the official IAU B1875.0 constellation boundaries.
    # The segment is a line of constant Right Ascension.
    ra_h = 1 + 39.0 / 60  # RA is 1h 39.0m
    p1_dec_deg = 0.0          # Declination of the first point is 0d 0'
    p2_dec_deg = 10 + 35.0 / 60 # Declination of the second point is 10d 35'

    # Step 2: Verify that the ecliptic crosses this boundary segment.
    # We use the obliquity of the ecliptic for the B1875.0 epoch.
    epsilon_deg = 23.4588 # Obliquity in degrees for B1875.0

    # Convert RA to radians for the calculation
    ra_deg = ra_h * 15
    ra_rad = math.radians(ra_deg)
    epsilon_rad = math.radians(epsilon_deg)

    # Calculate the declination of the ecliptic at the boundary's RA
    # The formula is tan(Dec) = tan(epsilon) * sin(RA)
    ecliptic_dec_rad = math.atan(math.tan(epsilon_rad) * math.sin(ra_rad))
    ecliptic_dec_deg = math.degrees(ecliptic_dec_rad)
    
    # Check if the ecliptic's declination is between the segment's endpoints
    is_crossing = p1_dec_deg < ecliptic_dec_deg < p2_dec_deg
    
    if not is_crossing:
        print("Error: The ecliptic does not cross the calculated boundary segment.")
        return

    # Step 3: Format the output as specified: XX YY ZZ, AA.BB
    # The "numerically lowest point" appears first. Since RA is the same,
    # we sort by declination.
    
    # Point 1 (lower declination)
    p1_ra_h, p1_ra_m, p1_ra_s = format_ra_hms(ra_h)
    p1_dec_str = format_dec_deg(p1_dec_deg).zfill(5)
    point1_str = f"{p1_ra_h} {p1_ra_m} {p1_ra_s}, {p1_dec_str}"

    # Point 2 (higher declination)
    p2_ra_h, p2_ra_m, p2_ra_s = format_ra_hms(ra_h)
    p2_dec_str = format_dec_deg(p2_dec_deg)
    point2_str = f"{p2_ra_h} {p2_ra_m} {p2_ra_s}, {p2_dec_str}"
    
    # The prompt requests "output each number in the final equation!".
    # This means showing the components of the final formatted string.
    # Final answer is the two points separated by a semicolon.
    final_answer_string = f"{point1_str}; {point2_str}"

    print(final_answer_string)


solve_boundary_crossing()
<<<01 39 00, 00.00; 01 39 00, 10.58>>>