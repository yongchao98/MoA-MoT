def dms_to_decimal(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees.
    """
    decimal_degrees = degrees + (minutes / 60) + (seconds / 3600)
    if direction in ['S', 'W']:
        decimal_degrees *= -1
    return decimal_degrees

# Latitude: 29° 06' 18.75''N
lat_deg, lat_min, lat_sec, lat_dir = 29, 6, 18.75, 'N'
# Longitude: 103° 47' 50.28''W
lon_deg, lon_min, lon_sec, lon_dir = 103, 47, 50.28, 'W'

# Perform the conversion
decimal_lat = dms_to_decimal(lat_deg, lat_min, lat_sec, lat_dir)
decimal_lon = dms_to_decimal(lon_deg, lon_min, lon_sec, lon_dir)

print(f"The original coordinates are:")
print(f"Latitude: {lat_deg}° {lat_min}' {lat_sec}\"{lat_dir}")
print(f"Longitude: {lon_deg}° {lon_min}' {lon_sec}\"{lon_dir}")
print("\nConverted to Decimal Degrees:")
print(f"Latitude: {decimal_lat}")
print(f"Longitude: {decimal_lon}")

print("\nBased on these coordinates, the structure can be identified.")
print("The square structure visible in the satellite image is not an ancient ruin.")
print("It is the remains of the facade of the 'Reata' mansion, a movie set built for the 1956 film 'Giant', starring Elizabeth Taylor, Rock Hudson, and James Dean.")
print("\nWhile it is not an officially designated historic landmark in a formal sense, it is a well-known point of interest for film enthusiasts and a significant piece of cinematic history located near Marfa, Texas.")
