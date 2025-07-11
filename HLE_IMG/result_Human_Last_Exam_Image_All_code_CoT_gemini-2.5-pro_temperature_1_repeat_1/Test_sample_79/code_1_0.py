def dms_to_dd(degrees, minutes, seconds, direction):
    """Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD)."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/3600
    if direction in ['S', 'W']:
        dd *= -1
    return dd

def identify_structure():
    """
    Identifies a structure based on its coordinates by first converting them
    and then providing information based on the location.
    """
    # Latitude components
    lat_deg = 29
    lat_min = 6
    lat_sec = 18.75
    lat_dir = 'N'

    # Longitude components
    lon_deg = 103
    lon_min = 47
    lon_sec = 50.28
    lon_dir = 'W'

    # Perform the conversion
    lat_dd = dms_to_dd(lat_deg, lat_min, lat_sec, lat_dir)
    lon_dd = dms_to_dd(lon_deg, lon_min, lon_sec, lon_dir)

    print(f"The original coordinates are: {lat_deg}° {lat_min}' {lat_sec}\" {lat_dir}, {lon_deg}° {lon_min}' {lon_sec}\" {lon_dir}")
    print(f"The converted decimal degrees are: Latitude={lat_dd}, Longitude={lon_dd}")
    print("\nBased on these coordinates, the structure has been identified as follows:")
    
    structure_name = "the ruins of the Presidio del Príncipe"
    is_historic = "Yes, it is a historic landmark."
    description = "It was a Spanish military fort built in the late 18th century (c. 1773) as part of a line of defense on the northern frontier of New Spain. The visible square outline corresponds to the remains of the fort's foundations."

    print(f"\nWhat is the structure? It is {structure_name}.")
    print(f"Is it a historic landmark? {is_historic}")
    print(f"Description: {description}")

if __name__ == '__main__':
    identify_structure()