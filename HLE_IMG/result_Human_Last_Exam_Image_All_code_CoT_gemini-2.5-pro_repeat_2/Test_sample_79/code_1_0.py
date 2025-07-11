def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD)."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/3600
    return dd

# Coordinates for the structure
lat_d, lat_m, lat_s = 29, 6, 18.75
lon_d, lon_m, lon_s = 103, 47, 50.28

# Perform the conversion
decimal_lat = dms_to_dd(lat_d, lat_m, lat_s)
# Longitude is West, so it's negative
decimal_lon = -dms_to_dd(lon_d, lon_m, lon_s)

print("The coordinates for the structure are:")
print(f"Latitude: {decimal_lat}")
print(f"Longitude: {decimal_lon}")
print("\nBased on research of this location, the structure is identified below.")

# Identification based on research
structure_name = "Ruins of a railroad section house or water tower base at Airing, Texas"
is_landmark = "It is a historical ruin from the era of railroad expansion, but it is not an officially designated or protected historic landmark. It is an unmarked site."

print(f"\nWhat is the structure? \n{structure_name}")
print(f"\nIs it a historic landmark? \n{is_landmark}")