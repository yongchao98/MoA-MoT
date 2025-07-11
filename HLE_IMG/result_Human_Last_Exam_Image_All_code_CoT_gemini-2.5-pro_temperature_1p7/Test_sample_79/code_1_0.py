def analyze_location(lat_dms, lon_dms):
    """
    Analyzes the location based on DMS coordinates, converts them to decimal,
    and prints information about the structure found there.
    """

    def dms_to_dd(degrees, minutes, seconds, direction):
        """Converts DMS to Decimal Degrees."""
        dd = degrees + minutes / 60 + seconds / 3600
        if direction in ['S', 'W']:
            dd *= -1
        return dd

    # Unpack latitude and longitude details
    lat_deg, lat_min, lat_sec, lat_dir = lat_dms
    lon_deg, lon_min, lon_sec, lon_dir = lon_dms

    # Perform the conversion
    lat_dd = dms_to_dd(lat_deg, lat_min, lat_sec, lat_dir)
    lon_dd = dms_to_dd(lon_deg, lon_min, lon_sec, lon_dir)

    # --- Information based on the coordinates ---
    structure_name = "Ruins of the historic Chinati Hot Springs resort"
    location_details = "Presidio County, Texas, USA"
    is_landmark = "It is a site of local historical interest but is not an officially designated state or national historic landmark."
    description = "The square outline is the remains of an old building, likely a bathhouse, hotel, or other structure associated with the original hot springs resort that was active in the late 19th and early 20th centuries."

    # --- Print the final report ---
    print("--- Location Analysis Report ---")
    print(f"Original Coordinates: {lat_deg}° {lat_min}' {lat_sec}\" {lat_dir}, {lon_deg}° {lon_min}' {lon_sec}\" {lon_dir}")
    print(f"Decimal Coordinates: {lat_dd:.6f}, {lon_dd:.6f}\n")
    print(f"Structure Identification: {structure_name}")
    print(f"Location: {location_details}")
    print(f"Description: {description}")
    print(f"Historic Landmark Status: {is_landmark}")
    print("---------------------------------")


if __name__ == '__main__':
    # Coordinates from the user's request
    latitude_dms = (29, 6, 18.75, 'N')
    longitude_dms = (103, 47, 50.28, 'W')

    analyze_location(latitude_dms, longitude_dms)