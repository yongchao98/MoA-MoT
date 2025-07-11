def dms_to_decimal(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees.
    """
    decimal_degrees = degrees + minutes/60 + seconds/3600
    if direction in ['S', 'W']:
        decimal_degrees *= -1
    return decimal_degrees

# Latitude: 29° 06' 18.75''N
lat_deg = 29
lat_min = 6
lat_sec = 18.75
lat_dir = 'N'

# Longitude: 103° 47' 50.28''W
lon_deg = 103
lon_min = 47
lon_sec = 50.28
lon_dir = 'W'

# Perform the conversion
decimal_lat = dms_to_decimal(lat_deg, lat_min, lat_sec, lat_dir)
decimal_lon = dms_to_decimal(lon_deg, lon_min, lon_sec, lon_dir)

# Construct Google Maps URL
google_maps_url = f"https://www.google.com/maps?q={decimal_lat},{decimal_lon}&z=17&t=h"

# Print the results
print("--- Location Identification ---")
print(f"Original Coordinates: {lat_deg}° {lat_min}' {lat_sec:.2f}''{lat_dir}, {lon_deg}° {lon_min}' {lon_sec:.2f}''{lon_dir}")
print(f"Converted to Decimal Degrees: {decimal_lat}, {decimal_lon}")
print("\n--- Structure Identification ---")
print("The structure at these coordinates is the remnant foundation of a building from the Marfa Army Airfield, a historic World War II-era airbase.")
print("The site is located near Marfa, Texas.")
print("\n--- Verification Link ---")
print("You can view the location on Google Maps (Satellite View):")
print(google_maps_url)
