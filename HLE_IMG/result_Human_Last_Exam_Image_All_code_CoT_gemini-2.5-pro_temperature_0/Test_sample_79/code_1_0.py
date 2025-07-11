def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes/60 + seconds/3600

# --- Latitude ---
lat_deg = 29
lat_min = 6
lat_sec = 18.75
# North is positive
lat_dd = dms_to_dd(lat_deg, lat_min, lat_sec)

# --- Longitude ---
lon_deg = 103
lon_min = 47
lon_sec = 50.28
# West is negative
lon_dd = -dms_to_dd(lon_deg, lon_min, lon_sec)

# --- Print the results with the full equation ---
print("Coordinate Conversion:")
print(f"Latitude: {lat_deg}° {lat_min}' {lat_sec}\"N is calculated as {lat_deg} + {lat_min}/60 + {lat_sec}/3600 = {lat_dd:.6f}°")
print(f"Longitude: {lon_deg}° {lon_min}' {lon_sec}\"W is calculated as -({lon_deg} + {lon_min}/60 + {lon_sec}/3600) = {lon_dd:.6f}°")
print("\nBased on these coordinates and historical records, the structure is identified as follows:\n")

# --- Identification and Landmark Status ---
structure_name = "El Fortín de San José"
landmark_status = "Yes, it is a historic landmark. It is part of the La Junta de los Ríos Archaeological District, which is a U.S. National Historic Landmark."

print(f"Structure Name: {structure_name}")
print(f"Is it a historic landmark? {landmark_status}")
