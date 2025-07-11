import math

def solve_constellation_boundary():
    """
    This function determines and prints the IAU boundary that the Vernal Equinox
    most recently crossed to enter the constellation Pisces.
    
    The logic is based on the fact that due to an intrusion of the constellation Cetus,
    the precessing Vernal Equinox first crossed from Aries into Cetus, and then later
    crossed from Cetus into Pisces. This latter boundary is the one we need.
    
    This boundary is a line of constant Right Ascension in the B1875.0 epoch.
    Data is from the official IAU 1930 "Délimitation scientifique des constellations".
    """
    
    # This boundary is located at Right Ascension 23h 46.2m for the B1875.0 epoch.
    ra_h = 23
    ra_m_decimal = 46.2
    
    # Convert RA from decimal minutes to minutes and seconds
    ra_m = int(ra_m_decimal)
    ra_s = int(round((ra_m_decimal - ra_m) * 60))

    # The boundary segment is defined by two declination points.
    # Point 1 (Southern endpoint): Declination -6 degrees 35 arcminutes
    dec1_deg = -6
    dec1_min = 35
    
    # Point 2 (Northern endpoint): Declination +5 degrees 41 arcminutes
    dec2_deg = 5
    dec2_min = 41
    
    # Convert declinations to decimal degrees, formatted to two decimal places
    dec1_decimal = dec1_deg - (dec1_min / 60.0)
    dec2_decimal = dec2_deg + (dec2_min / 60.0)
    
    # The "numerically lowest" point is the one with the lower declination value.
    # Point_lower corresponds to dec1.
    # Point_upper corresponds to dec2.

    # Print the explanation and the final result
    print("The boundary marker where the Vernal Equinox most recently crossed into Pisces is a line segment separating Pisces from Cetus.")
    print("The two endpoints defining this line segment are calculated and formatted below.")
    print("Lower Point (RA, Dec):")
    print(f"RA: {ra_h}h {ra_m}m {ra_s}s, Declination: {dec1_decimal:.2f} degrees")
    print("\nUpper Point (RA, Dec):")
    print(f"RA: {ra_h}h {ra_m}m {ra_s}s, Declination: {dec2_decimal:.2f} degrees")
    print("\nFinal formatted answer:")
    
    # Format the final output string as "XX YY ZZ, AA.BB; XX YY ZZ, AA.BB"
    # The point with the lower declination (-6.58) comes first.
    final_output = (f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec1_decimal:+.2f}; "
                    f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec2_decimal:+.2f}")
    
    print(final_output)

solve_constellation_boundary()