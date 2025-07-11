def identify_structure_from_coords():
    """
    Converts coordinates from DMS to DD and identifies the structure at that location.
    """
    # Latitude: 29° 06' 18.75''N
    lat_deg = 29
    lat_min = 6
    lat_sec = 18.75

    # Longitude: 103° 47' 50.28''W
    lon_deg = 103
    lon_min = 47
    lon_sec = 50.28

    print("--- Coordinate Conversion ---")

    # --- Latitude Calculation ---
    print("\nCalculating Latitude (N):")
    # Equation: DD = Degrees + (Minutes / 60) + (Seconds / 3600)
    lat_dd = lat_deg + (lat_min / 60) + (lat_sec / 3600)
    print(f"Equation: {lat_deg} + ({lat_min} / 60) + ({lat_sec} / 3600)")
    print(f"Result: {lat_deg} + {lat_min / 60:.6f} + {lat_sec / 3600:.6f} = {lat_dd:.6f}")
    
    # --- Longitude Calculation ---
    print("\nCalculating Longitude (W):")
    # Equation: DD = -(Degrees + (Minutes / 60) + (Seconds / 3600))
    lon_dd = -1 * (lon_deg + (lon_min / 60) + (lon_sec / 3600))
    print(f"Equation: -({lon_deg} + ({lon_min} / 60) + ({lon_sec} / 3600))")
    print(f"Result: -({lon_deg} + {lon_min / 60:.6f} + {lon_sec / 3600:.6f}) = {lon_dd:.6f}")
    
    print("\n--- Identification ---")
    print("The structure at these coordinates is the ruins of the Presidio de San Vicente.")
    print("\nIs it a historic landmark?")
    print("Yes, it is a significant historic landmark. It was a Spanish fort (presidio) built around 1774 as part of a line of defenses along the northern frontier of New Spain to protect against Apache raids. It is located in Chihuahua, Mexico, near the Rio Grande.")

if __name__ == '__main__':
    identify_structure_from_coords()