def dms_to_dd(name, degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD).
    Also prints the calculation process.
    """
    decimal_degrees = degrees + (minutes / 60) + (seconds / 3600)
    
    # West and South coordinates are negative
    if direction in ['W', 'S']:
        decimal_degrees *= -1
    
    print(f"Calculating Decimal Degrees for {name}:")
    # Building the equation string for output
    equation = ""
    if direction in ['W', 'S']:
        equation += "-("
        
    equation += f"{degrees}° + {minutes}'/60 + {seconds}\"/3600"
    
    if direction in ['W', 'S']:
        equation += ")"
    
    equation += f" = {decimal_degrees:.6f}"
    
    print(equation)
    return decimal_degrees

# Given coordinates: 29° 3' 28.15''N and 103° 48' 11.84''W
lat_d, lat_m, lat_s = 29, 3, 28.15
lon_d, lon_m, lon_s = 103, 48, 11.84

# Perform the conversions
print("--- Coordinate Conversion ---")
decimal_lat = dms_to_dd("Latitude", lat_d, lat_m, lat_s, 'N')
decimal_lon = dms_to_dd("Longitude", lon_d, lon_m, lon_s, 'W')

print("\n--- Final Coordinates ---")
print(f"The decimal coordinates for the rock carving are approximately:")
print(f"Latitude: {decimal_lat:.6f}")
print(f"Longitude: {decimal_lon:.6f}")
print("\nUsing these coordinates, you can verify the location is about 10 miles from the Rio Bravo near Lajitas and about 20 miles from the Chisos Mountains.")
