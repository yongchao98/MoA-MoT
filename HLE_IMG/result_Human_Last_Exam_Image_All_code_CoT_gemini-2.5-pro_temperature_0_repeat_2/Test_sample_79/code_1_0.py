def identify_structure_from_coords():
    """
    Identifies a structure based on given geographical coordinates and
    provides details about it.
    """
    # --- Input Coordinates in Degrees, Minutes, Seconds ---
    lat_d, lat_m, lat_s = 29, 6, 18.75
    lon_d, lon_m, lon_s = 103, 47, 50.28

    # --- Step 1: Convert DMS to Decimal Degrees ---
    # Latitude is North (positive)
    lat_dd = lat_d + lat_m / 60 + lat_s / 3600
    # Longitude is West (negative)
    lon_dd = -1 * (lon_d + lon_m / 60 + lon_s / 3600)

    # --- Step 2: Identify the structure and its details ---
    # Based on the coordinates, the structure is identified as Prada Marfa.
    structure_name = "Prada Marfa"
    structure_type = "Permanent art installation"
    year_built = 2005
    is_historic = False

    # --- Step 3: Print the results ---
    print("--- Analysis of the Structure at the Given Coordinates ---")
    print(f"\nCoordinates (DMS): {lat_d}° {lat_m}' {lat_s}\" N, {lon_d}° {lon_m}' {lon_s}\" W")

    print("\nCoordinate Conversion to Decimal Degrees:")
    # The prompt requires showing the numbers in the final equation
    print(f"Latitude Equation: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {lat_dd:.6f}")
    print(f"Longitude Equation: -1 * ({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600)) = {lon_dd:.6f}")

    print("\n--- Identification ---")
    print(f"The square structure in the image is not an ancient building.")
    print(f"It is a well-known art installation called '{structure_name}'.")

    print("\n--- Details ---")
    print(f"  - What it is: {structure_type}")
    print(f"  - Year of creation: {year_built}")
    print(f"  - Is it a historic landmark? No, it is a contemporary cultural landmark, not a historic one.")
    print("    It was created in 2005 by artists Elmgreen & Dragset. While it's a famous landmark, it's too recent to be officially 'historic'.")

if __name__ == '__main__':
    identify_structure_from_coords()