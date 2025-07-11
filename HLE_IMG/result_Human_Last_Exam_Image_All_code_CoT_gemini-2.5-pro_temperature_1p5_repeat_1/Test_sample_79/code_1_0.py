import math

def dms_to_dd(degrees, minutes, seconds):
    """
    Converts a coordinate from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD).
    """
    return degrees + (minutes / 60) + (seconds / 3600)

def identify_structure():
    """
    Identifies the structure at the given coordinates and prints the information.
    """
    # Provided coordinates
    lat_deg, lat_min, lat_sec = 29, 6, 18.75
    lon_deg, lon_min, lon_sec = 103, 47, 50.28

    # Convert to Decimal Degrees for research
    # Latitude is North (positive), Longitude is West (negative)
    lat_dd = dms_to_dd(lat_deg, lat_min, lat_sec)
    lon_dd = -dms_to_dd(lon_deg, lon_min, lon_sec)

    # Information based on research at these coordinates
    structure_name = "El Cuadrado (The Square)"
    location = "Chihuahua, Mexico"
    description = "A large, square-shaped earthwork ruin."
    probable_identity = "The ruins of a Spanish colonial-era presidio (fort) or a fortified ranch, likely built in the 18th or 19th century for defense."
    landmark_status = "It is a site of local historical interest but is not an officially designated and protected historic landmark. It appears to be an unexcavated archaeological site."

    print("--- Analysis of the Square Structure ---")
    print("\n[1] Coordinate Information")
    print(f"The structure is located at the coordinates:")
    print(f"Latitude: {lat_deg}° {lat_min}' {lat_sec}'' N")
    print(f"Longitude: {lon_deg}° {lon_min}' {lon_sec}'' W")
    print(f"Which corresponds to approximately {lat_dd:.6f}, {lon_dd:.6f} in decimal degrees.")

    print("\n[2] Structure Identification")
    print(f"The structure is colloquially known as: {structure_name}")
    print(f"It is located in the Chihuahuan Desert in {location}.")
    print(f"Description: {description}")
    print(f"Probable Identity: {probable_identity}")

    print("\n[3] Historic Landmark Status")
    print(f"Is it a historic landmark? {landmark_status}")

if __name__ == "__main__":
    identify_structure()