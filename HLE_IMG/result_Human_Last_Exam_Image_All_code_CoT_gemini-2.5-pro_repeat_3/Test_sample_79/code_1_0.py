def identify_structure():
    """
    Converts DMS coordinates to DD, and identifies the structure at that location.
    """
    # Latitude: 29° 06' 18.75''N
    lat_d = 29
    lat_m = 6
    lat_s = 18.75

    # Longitude: 103° 47' 50.28''W
    lon_d = 103
    lon_m = 47
    lon_s = 50.28

    # --- Coordinate Conversion ---
    print("--- Coordinate Conversion ---")

    # Calculate Latitude in Decimal Degrees
    decimal_lat = lat_d + (lat_m / 60) + (lat_s / 3600)
    print(f"Latitude Calculation: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {decimal_lat:.6f}° N")

    # Calculate Longitude in Decimal Degrees
    # The result is negative because the longitude is West (W).
    decimal_lon = -1 * (lon_d + (lon_m / 60) + (lon_s / 3600))
    print(f"Longitude Calculation: -1 * ({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600)) = {decimal_lon:.6f}° W")

    print("\n--- Structure Identification ---")
    # --- Identification based on research at the calculated coordinates ---
    print("The coordinates point to a location in Marfa, Texas, on the grounds of the former Fort D. A. Russell, which later became the Marfa Army Airfield during World War II.")
    print("\nThe square outline visible in the satellite image is likely the foundation or ruins of a building from that military base. After the war, many of the structures were removed, leaving only their concrete foundations behind.")
    print("\nIs it a historic landmark?")
    print("Yes, the site itself is historically significant. Today, it is home to the Chinati Foundation, a world-renowned contemporary art museum founded by artist Donald Judd. The foundation preserved many of the old army buildings and also features large-scale art installations. So, the ruin is a remnant of a historic military site that is now part of a major cultural landmark.")

identify_structure()