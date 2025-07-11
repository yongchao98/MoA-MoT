# --- Coordinate Conversion ---

# Latitude: 29° 06' 18.75''N
lat_deg = 29
lat_min = 6
lat_sec = 18.75

# Longitude: 103° 47' 50.28''W
lon_deg = 103
lon_min = 47
lon_sec = 50.28

# Calculate Decimal Degrees
# The sign for North latitude is positive
decimal_lat = lat_deg + (lat_min / 60) + (lat_sec / 3600)
# The sign for West longitude is negative
decimal_lon = -1 * (lon_deg + (lon_min / 60) + (lon_sec / 3600))

# Print the explanation and results
print("Converting the coordinates from DMS (Degrees, Minutes, Seconds) to DD (Decimal Degrees).")
print("\n--- Latitude Conversion (29° 06' 18.75''N) ---")
print(f"Formula: Degrees + (Minutes / 60) + (Seconds / 3600)")
print(f"Calculation: {lat_deg} + ({lat_min} / 60) + ({lat_sec} / 3600) = {decimal_lat:.6f}")


print("\n--- Longitude Conversion (103° 47' 50.28''W) ---")
print(f"Formula: -1 * (Degrees + (Minutes / 60) + (Seconds / 3600))")
print(f"Calculation: -1 * ({lon_deg} + ({lon_min} / 60) + ({lon_sec} / 3600)) = {decimal_lon:.6f}")

print(f"\nThe coordinates point to approximately {decimal_lat:.6f}, {decimal_lon:.6f}.")

print("\n--- Identification ---")
print("This structure is the ruin of the Presidio de San Vicente.")
print("\n--- Historical Significance ---")
print("Yes, it is a significant historic landmark. It was a Spanish fort built in the 1770s to defend the northern frontier of New Spain. It is listed on the National Register of Historic Places.")