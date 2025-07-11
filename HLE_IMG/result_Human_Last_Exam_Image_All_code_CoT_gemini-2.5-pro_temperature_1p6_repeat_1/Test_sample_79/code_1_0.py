def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + (minutes / 60) + (seconds / 3600)

def main():
    """
    Identifies a structure from given coordinates and provides historical context.
    """
    # Coordinates from the user query
    lat_d, lat_m, lat_s = 29, 6, 18.75
    lon_d, lon_m, lon_s = 103, 47, 50.28

    # Convert to Decimal Degrees
    decimal_lat = dms_to_dd(lat_d, lat_m, lat_s)
    # Longitude is West, so it's negative
    decimal_lon = -dms_to_dd(lon_d, lon_m, lon_s)

    # Based on research of the coordinates (29.1052 N, -103.7973 W)
    structure_name = "Presidio de San Vicente"
    is_landmark = "Yes"
    location = "Coahuila, Mexico, near the Rio Grande"
    description = (
        "This square structure consists of the ruins of the Presidio de San Vicente. "
        "It was a Spanish colonial military fort constructed around 1773 as part of a line of "
        "defenses along the northern frontier of New Spain to protect against Apache raids."
    )
    historical_significance = (
        "It is a significant historic landmark and an important archaeological site. "
        "The presidio represents the Spanish Empire's efforts to control and settle the "
        "region in the late 18th century. The site is under the protection of Mexico's "
        "National Institute of Anthropology and History (INAH)."
    )

    print("--- Identification of the Structure ---")
    print(f"Coordinates (DMS): {lat_d}° {lat_m}' {lat_s}\"N, {lon_d}° {lon_m}' {lon_s}\"W")
    print(f"Coordinates (Decimal): {decimal_lat:.4f}, {decimal_lon:.4f}\n")

    print(f"What is this structure? \n> {structure_name}\n")
    print(f"Location: \n> {location}\n")
    print(f"Description: \n> {description}\n")

    print(f"Is it a historic landmark? \n> {is_landmark}\n")
    print(f"Historical Significance: \n> {historical_significance}")

if __name__ == "__main__":
    main()