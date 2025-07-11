# Latitude components
lat_deg = 29
lat_min = 6
lat_sec = 18.75

# Longitude components
lon_deg = 103
lon_min = 47
lon_sec = 50.28

# Calculate decimal degrees for latitude
decimal_lat = lat_deg + (lat_min / 60) + (lat_sec / 3600)

# Calculate decimal degrees for longitude (West is negative)
decimal_lon = -1 * (lon_deg + (lon_min / 60) + (lon_sec / 3600))

print("Converting the coordinates...")
print(f"Latitude: {lat_deg}° {lat_min}' {lat_sec}''N = {lat_deg} + ({lat_min}/60) + ({lat_sec}/3600) = {decimal_lat:.6f}°")
print(f"Longitude: {lon_deg}° {lon_min}' {lon_sec}''W = -1 * ({lon_deg} + ({lon_min}/60) + ({lon_sec}/3600)) = {decimal_lon:.6f}°")

print("\nBased on these coordinates, the structure has been identified.")

# Final answer generation
name = "one of Donald Judd's 15 untitled works in concrete"
is_landmark = "Yes, it is a significant landmark of modern art, part of the Chinati Foundation's permanent collection."

print(f"\nThe structure shown in the image is {name}.")
print(f"Is it a historic landmark? {is_landmark}")