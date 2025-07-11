def solve_boundary_question():
    """
    Calculates and prints the formatted boundary coordinates between Pisces and Aries.
    """
    # The boundary between Aries and Pisces that crosses the celestial equator is a
    # line of constant Right Ascension. The coordinates for the two vertices defining
    # this segment are based on the IAU B1875.0 epoch definitions.

    # Right Ascension is constant for both points
    RA_H = 23
    RA_M = 59
    RA_S = 19.3

    # Vertex 1 (Northern point)
    DEC1_D = 5
    DEC1_M = 28
    DEC1_S = 22

    # Vertex 2 (Southern point)
    DEC2_D = -2
    DEC2_M = 34
    DEC2_S = 13

    # Convert the declination for the northern point to decimal degrees.
    # Dec = D + M/60 + S/3600
    dec1_decimal = DEC1_D + (DEC1_M / 60.0) + (DEC1_S / 3600.0)

    # Convert the declination for the southern point to decimal degrees.
    # For negative declinations given as -D M S, the calculation is -(D + M/60 + S/3600)
    dec2_decimal = DEC2_D - (DEC2_M / 60.0) - (DEC2_S / 3600.0)

    # The problem requires the numerically lowest point (southernmost) to be listed first.
    # The output format is "XX YY ZZ, AA.BB; XX YY ZZ, AA.BB"
    # We will print the final string assembled from its constituent parts.
    
    print(
        f"{RA_H} {RA_M} {RA_S}, {dec2_decimal:.2f}; "
        f"{RA_H} {RA_M} {RA_S}, {dec1_decimal:.2f}"
    )

solve_boundary_question()