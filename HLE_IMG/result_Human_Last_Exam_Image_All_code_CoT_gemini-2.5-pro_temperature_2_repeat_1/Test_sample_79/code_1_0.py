def dms_to_dd(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD).
    """
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

# --- Latitude ---
lat_deg = 29
lat_min = 6
lat_sec = 18.75
lat_dir = 'N'

# --- Longitude ---
lon_deg = 103
lon_min = 47
lon_sec = 50.28
lon_dir = 'W'

# --- Conversion ---
decimal_lat = dms_to_dd(lat_deg, lat_min, lat_sec, lat_dir)
decimal_lon = dms_to_dd(lon_deg, lon_min, lon_sec, lon_dir)

print(f"Original Latitude: {lat_deg}° {lat_min}' {lat_sec}\" {lat_dir}")
print(f"Original Longitude: {lon_deg}° {lon_min}' {lon_sec}\" {lon_dir}")
print("-" * 30)
print(f"Converted to Decimal Degrees:")
print(f"Latitude: {decimal_lat}")
print(f"Longitude: {decimal_lon}")

print("\n--- Structure Identification ---")
print("The structure at these coordinates is a large, modern earthwork, often referred to as 'The Square' or the 'Cíbolo Creek Square'.")
print("Despite its appearance, it is not an ancient structure. It is believed to be a contemporary land art installation, a common feature in this region famous for the minimalist art scene of Marfa, Texas.")
print("\nIs it a historic landmark?")
print("No, this specific square structure is not an officially recognized historic landmark. A true historic site, the ruins of the 18th-century Spanish 'Presidio de San Carlos', is located nearby but is a separate location.")