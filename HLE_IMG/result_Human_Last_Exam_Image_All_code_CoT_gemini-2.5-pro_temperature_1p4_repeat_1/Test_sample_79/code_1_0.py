def dms_to_dd(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD).
    """
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

# Latitude: 29° 06' 18.75''N
lat_d = 29
lat_m = 6
lat_s = 18.75
lat_dir = 'N'

# Longitude: 103° 47' 50.28''W
lon_d = 103
lon_m = 47
lon_s = 50.28
lon_dir = 'W'

# Perform the conversion
decimal_latitude = dms_to_dd(lat_d, lat_m, lat_s, lat_dir)
decimal_longitude = dms_to_dd(lon_d, lon_m, lon_s, lon_dir)

# Print the results and the identification
print(f"Identifying the structure at the given coordinates:")
print(f"Original Latitude (DMS): {lat_d}° {lat_m}' {lat_s}\"{lat_dir}")
print(f"Original Longitude (DMS): {lon_d}° {lon_m}' {lon_s}\"{lon_dir}")
print("-" * 30)
print(f"Converted Latitude (DD): {decimal_latitude}")
print(f"Converted Longitude (DD): {decimal_longitude}")
print("-" * 30)
print("Identification: The structure is the ruins of Camp Cienega Springs.")
print("Historic Landmark: Yes, it is a historic site from the late 19th-century U.S. Army presence in West Texas.")
