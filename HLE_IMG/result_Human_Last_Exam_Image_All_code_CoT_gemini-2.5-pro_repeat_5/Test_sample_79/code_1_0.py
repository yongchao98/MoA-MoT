def dms_to_dd(degrees, minutes, seconds, direction):
    """Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD)."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

# Given coordinates
lat_d, lat_m, lat_s, lat_dir = 29, 6, 18.75, 'N'
lon_d, lon_m, lon_s, lon_dir = 103, 47, 50.28, 'W'

# Perform the conversion
lat_dd = dms_to_dd(lat_d, lat_m, lat_s, lat_dir)
lon_dd = dms_to_dd(lon_d, lon_m, lon_s, lon_dir)

print("Step 1: Convert the coordinates from DMS to Decimal Degrees.")
print("The formula for latitude (N) is: DD = Degrees + (Minutes / 60) + (Seconds / 3600)")
print(f"Latitude = {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {lat_dd:.6f}")
print("\nThe formula for longitude (W) is: DD = -(Degrees + (Minutes / 60) + (Seconds / 3600))")
print(f"Longitude = -({lon_d} + ({lon_m} / 60) + ({lon_s} / 3600)) = {lon_dd:.6f}")

print("\nStep 2: Research the location using the decimal coordinates.")
print("Based on research of the coordinates and the satellite image, the identity of the structure can be determined.")

print("\n--- Identification ---")
print("What is the structure? The square outline corresponds to the ruins of the Presidio de San Vicente.")
print("This was a Spanish fort built around 1774 as part of a defensive line along the northern frontier of New Spain.")

print("\nIs it a historic landmark? Yes, it is a historic landmark.")
print("The site is located within Big Bend Ranch State Park in Texas and is listed on the National Register of Historic Places.")