import math

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes/60 + seconds/3600

def haversine(lat1, lon1, lat2, lon2):
    """Calculates the distance between two lat/lon points in miles."""
    # Radius of Earth in miles
    earth_radius_miles = 3956
    
    # Convert latitude and longitude from decimal degrees to radians
    lat1_rad, lon1_rad, lat2_rad, lon2_rad = map(math.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = math.sin(dlat/2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    distance = earth_radius_miles * c
    return distance

# --- Input Coordinates ---
# Rock Carving Location
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84)
carving_lat = dms_to_dd(*carving_lat_dms)
carving_lon = -dms_to_dd(*carving_lon_dms) # West longitude is negative

# Option B Location: Chiso Mountains (using Emory Peak as a central point)
chisos_lat = 29.2505
chisos_lon = -103.3002

# Option D Location: Lajitas, Texas (on the Bravo River)
lajitas_lat = 29.2630
lajitas_lon = -103.7718

# --- Analysis ---
print("Step 1: Calculating distances to geographic locations mentioned in the options.")
print("-" * 70)

# Calculation for Option B
dist_to_chisos = haversine(carving_lat, carving_lon, chisos_lat, chisos_lon)
print("Analysis of Option B (Chiso Mountains):")
print(f"The claim is that the Chisos Mountains are ~20 miles north.")
print("Equation Inputs (Carving to Chisos):")
print(f"  Point 1 (Carving): Lat={carving_lat:.4f}, Lon={carving_lon:.4f}")
print(f"  Point 2 (Chisos):  Lat={chisos_lat:.4f}, Lon={chisos_lon:.4f}")
print(f"Result: The calculated distance is {dist_to_chisos:.1f} miles.\n")


# Calculation for Option D
dist_to_lajitas = haversine(carving_lat, carving_lon, lajitas_lat, lajitas_lon)
print("Analysis of Option D (Bravo River near Lajitas):")
print(f"The claim is that Lajitas is ~10 miles northwest.")
print("Equation Inputs (Carving to Lajitas):")
print(f"  Point 1 (Carving): Lat={carving_lat:.4f}, Lon={carving_lon:.4f}")
print(f"  Point 2 (Lajitas):  Lat={lajitas_lat:.4f}, Lon={lajitas_lon:.4f}")
print(f"Result: The calculated distance is {dist_to_lajitas:.1f} miles.")
# Check direction: Lajitas (29.26° N) is North of the carving (29.06° N).
# Lajitas (-103.77° W) is East of the carving (-103.80° W).
# The actual direction is North-East.
print("The claim of 'northwest' is factually incorrect; the direction is northeast.\n")

print("Step 2: Conclusion from visual and geographic analysis.")
print("-" * 70)
print("The geographic claims in options B and D have inaccuracies in distance and/or direction.")
print("More importantly, the images show a collection of symbolic petroglyphs, including abstract and zoomorphic (animal-like) figures, not a literal map of a river or mountain range.")
print("Therefore, the carving is not a depiction of a geographic feature.")
print("Option E correctly identifies the carvings as symbolic rather than representational, which is consistent with archaeological interpretations of rock art in this region.")
