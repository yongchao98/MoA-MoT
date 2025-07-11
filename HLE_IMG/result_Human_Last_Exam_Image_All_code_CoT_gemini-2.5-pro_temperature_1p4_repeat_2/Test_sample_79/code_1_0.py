def identify_structure():
    """
    This script converts the given coordinates and identifies the structure.
    It then prints the details about the landmark.
    """

    # Given coordinates
    lat_d, lat_m, lat_s = 29, 6, 18.75
    lon_d, lon_m, lon_s = 103, 47, 50.28

    # Perform and print the conversion calculations
    lat_decimal = lat_d + (lat_m / 60) + (lat_s / 3600)
    lon_decimal = -1 * (lon_d + (lon_m / 60) + (lon_s / 3600))

    print("--- Coordinate Conversion ---")
    print(f"Latitude calculation: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {lat_decimal:.6f}")
    print(f"Longitude calculation: -1 * ({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600)) = {lon_decimal:.6f}")
    print("-" * 29)
    print("\n--- Structure Identification ---")
    
    # Based on the coordinates, the structure is identified.
    structure_name = "Presidio de San Vicente"
    is_historic_landmark = True
    description = (
        "The square structure shown in the image at the provided coordinates is the ruins of the "
        "Presidio de San Vicente. It is located in what is now Big Bend National Park, Texas.\n\n"
        "It was a Spanish fort, or 'presidio', established around 1774 as part of a line of "
        "defenses along the northern frontier of New Spain to protect against raids by groups "
        "like the Apache. The visible square outline corresponds to the remnants of the fort's "
        "adobe walls."
    )

    print(f"What is this structure? It is the {structure_name}.")
    print(f"\nIs it a historic landmark? {'Yes, it is a significant historic and archaeological site.' if is_historic_landmark else 'No.'}")
    print(f"\nDescription:\n{description}")

if __name__ == '__main__':
    identify_structure()
