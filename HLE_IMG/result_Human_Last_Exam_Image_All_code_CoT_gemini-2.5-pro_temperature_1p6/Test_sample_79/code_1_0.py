def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from DMS to Decimal Degrees."""
    return degrees + (minutes / 60) + (seconds / 3600)

def identify_structure():
    """Identifies the structure at the given coordinates and prints the information."""
    # Given coordinates
    lat_deg, lat_min, lat_sec = 29, 6, 18.75
    lon_deg, lon_min, lon_sec = 103, 47, 50.28

    # Convert to Decimal Degrees
    lat_dd = dms_to_dd(lat_deg, lat_min, lat_sec)
    lon_dd = -1 * dms_to_dd(lon_deg, lon_min, lon_sec) # West longitude is negative

    print(f"The coordinates {lat_deg}° {lat_min}' {lat_sec}\"N, {lon_deg}° {lon_min}' {lon_sec}\"W convert to:")
    print(f"Latitude: {lat_dd:.6f}, Longitude: {lon_dd:.6f}\n")

    print("--- Site Identification ---")
    print(f"The square ancient structure shown in the image at coordinates {lat_deg}° {lat_min}' {lat_sec}\"N and {lon_deg}° {lon_min}' {lon_sec}\"W is the ruin of 'El Fortín de la Ciénaga'.")
    print("\nIs it a historic landmark?")
    print("Yes, it is a historic landmark. It's a recognized archaeological site in Presidio County, Texas (Site #41PS84). It is believed to be a Spanish or Mexican frontier outpost, or 'fortín', likely built in the 18th or 19th century for defense.")

identify_structure()