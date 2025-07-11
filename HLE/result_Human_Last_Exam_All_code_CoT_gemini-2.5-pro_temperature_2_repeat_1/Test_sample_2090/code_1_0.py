import math

def solve():
    """
    This function finds and formats the IAU constellation boundary segment
    that the Vernal Equinox crossed when it most recently entered Pisces.
    """
    
    # Based on astronomical data (IAU 1930 boundaries, B1875.0 epoch), the
    # boundary segment is a line of constant Right Ascension that crosses the
    # celestial equator (Dec=0) on the eastern side of Pisces.
    
    # The vertices of this line segment are:
    p1_ra_h, p1_ra_m, p1_ra_s = 1, 42, 0
    p1_dec_d, p1_dec_m = -1, 5
    
    p2_ra_h, p2_ra_m, p2_ra_s = 1, 42, 0
    p2_dec_d, p2_dec_m = 1, 5
    
    # Calculate declination in decimal degrees.
    # The southern point's declination is negative.
    dec1 = p1_dec_d - (p1_dec_m / 60.0)
    # The northern point's declination is positive.
    dec2 = p2_dec_d + (p2_dec_m / 60.0)
    
    # The problem asks to output each number in the final result.
    # Let's print the components clearly.
    print("The boundary is defined by a line segment connecting two points.")
    print("\nPoint 1 (Southern vertex):")
    print(f"Right Ascension: {p1_ra_h}h {p1_ra_m}m {p1_ra_s}s")
    print(f"Declination: {p1_dec_d} degrees {abs(p1_dec_m)} minutes")
    
    print("\nPoint 2 (Northern vertex):")
    print(f"Right Ascension: {p2_ra_h}h {p2_ra_m}m {p2_ra_s}s")
    print(f"Declination: {p2_dec_d} degrees {p2_dec_m} minutes")
    
    # Format the points into the specified string format: "XX YY ZZ, AA.BB"
    # The points are ordered by declination (lowest first).
    p1_formatted = f"{p1_ra_h:02d} {p1_ra_m:02d} {p1_ra_s:02d}, {dec1:.2f}"
    p2_formatted = f"{p2_ra_h:02d} {p2_ra_m:02d} {p2_ra_s:02d}, {dec2:.2f}"
    
    # Combine the two points with a semicolon.
    final_answer = f"{p1_formatted}; {p2_formatted}"
    
    print("\nFinal formatted answer string:")
    print(final_answer)
    
    # Output the final answer in the required format for the system.
    print(f"<<<{final_answer}>>>")

solve()