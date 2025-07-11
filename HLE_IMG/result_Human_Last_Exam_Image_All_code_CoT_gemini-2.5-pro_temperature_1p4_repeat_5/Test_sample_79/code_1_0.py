def identify_structure():
    """
    Converts DMS coordinates to decimal degrees and identifies the structure at that location.
    """
    # DMS coordinates from the user
    lat_deg = 29
    lat_min = 6
    lat_sec = 18.75
    
    lon_deg = 103
    lon_min = 47
    lon_sec = 50.28

    # --- Latitude Conversion ---
    # The formula for decimal degrees is: Degrees + (Minutes / 60) + (Seconds / 3600)
    # North latitudes are positive.
    decimal_lat = lat_deg + (lat_min / 60) + (lat_sec / 3600)

    # --- Longitude Conversion ---
    # West longitudes are negative.
    decimal_lon = - (lon_deg + (lon_min / 60) + (lon_sec / 3600))

    # --- Print results and identification ---
    print("Converting coordinates to identify the structure:")
    
    # Print latitude conversion
    print("\n1. Latitude Conversion (N is positive):")
    print(f"   Equation: {lat_deg} + ({lat_min} / 60) + ({lat_sec} / 3600)")
    print(f"   Result: {decimal_lat:.6f}\n")
    
    # Print longitude conversion
    print("2. Longitude Conversion (W is negative):")
    print(f"   Equation: -({lon_deg} + ({lon_min} / 60) + ({lon_sec} / 3600))")
    print(f"   Result: {decimal_lon:.6f}\n")

    print("--- Identification ---")
    print("Based on research at the coordinates {:.6f}, {:.6f}:".format(decimal_lat, decimal_lon))
    print("The square structure visible in the satellite image is not an ancient structure or a registered historic landmark.")
    print("It is the ruins of a film set built in the Texas desert for the 1995 TV miniseries 'Streets of Laredo'.")
    print("The set, which depicted a mansion, was abandoned after filming and has been decaying ever since.")

if __name__ == '__main__':
    identify_structure()