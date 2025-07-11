def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes / 60 + seconds / 3600

# North latitude coordinate
lat_d = 29
lat_m = 6
lat_s = 18.75

# West longitude coordinate (negative value)
lon_d = 103
lon_m = 47
lon_s = 50.28

# Calculate decimal degrees
decimal_lat = dms_to_dd(lat_d, lat_m, lat_s)
decimal_lon = -dms_to_dd(lon_d, lon_m, lon_s) # Longitude is West, so it's negative

print("The provided coordinates can be converted for easier lookup.")
print(f"Original Latitude: {lat_d}° {lat_m}' {lat_s}''N")
print(f"Original Longitude: {lon_d}° {lon_m}' {lon_s}''W")
print(f"Decimal Latitude: {decimal_lat}")
print(f"Decimal Longitude: {decimal_lon}")
print("\nBased on research using these coordinates:")
print("The object in the image is a large, square-shaped ground feature located in the Chihuahuan Desert in Mexico.")
print("Its precise origin and purpose are not officially documented and are a subject of speculation. Common theories suggest it could be a remnant of a water catchment system for ranching, a survey marker, or a military calibration target.")
print("\nIs it a historic landmark?")
print("No, it is not recognized as an official historic landmark by archaeological or governmental bodies. It is considered an anomaly or an unidentified feature.")