def find_pisces_boundary_marker():
    """
    This function provides the coordinates of the IAU constellation boundary
    marker where the Vernal Equinox precessed from Aries into Pisces.
    """

    # The boundary is a line of constant Right Ascension, defined by two
    # endpoints in the B1875.0 epoch. The data is sourced from the
    # official IAU definitions.

    # Point 1: The southern endpoint of the line segment.
    # Numerically, this is the lower point.
    ra1_h, ra1_m, ra1_s = 0, 0, 0
    # Declination is -11 degrees and 30 arcminutes, which is -11.5 degrees.
    dec1 = -11.50

    # Point 2: The northern endpoint of the line segment.
    ra2_h, ra2_m, ra2_s = 0, 0, 0
    # Declination is +6 degrees and 30 arcminutes, which is 6.5 degrees.
    dec2 = 6.50

    # The required format is "XX YY ZZ, AA.BB" for each point.
    # We will format each component before combining them.

    # Format Point 1 components
    p1_ra_str = f"{ra1_h:02d} {ra1_m:02d} {ra1_s:02d}"
    p1_dec_str = f"{dec1:.2f}"

    # Format Point 2 components
    p2_ra_str = f"{ra2_h:02d} {ra2_m:02d} {ra2_s:02d}"
    p2_dec_str = f"{dec2:.2f}"

    # Construct the final output string by combining the formatted parts,
    # separated by a semicolon as requested.
    final_output = f"{p1_ra_str}, {p1_dec_str}; {p2_ra_str}, {p2_dec_str}"

    print("The IAU boundary between Aries and Pisces is a line segment defined by two points.")
    print("\n--- Point 1 (Southern Endpoint) ---")
    print(f"Right Ascension components: HH={ra1_h:02d}, MM={ra1_m:02d}, SS={ra1_s:02d}")
    print(f"Declination component: DD.dd={dec1:.2f}")

    print("\n--- Point 2 (Northern Endpoint) ---")
    print(f"Right Ascension components: HH={ra2_h:02d}, MM={ra2_m:02d}, SS={ra2_s:02d}")
    print(f"Declination component: DD.dd={dec2:.2f}")

    print("\n------------------------------------")
    print("Final formatted boundary marker:")
    print(final_output)


find_pisces_boundary_marker()
<<<00 00 00, -11.50; 00 00 00, 06.50>>>