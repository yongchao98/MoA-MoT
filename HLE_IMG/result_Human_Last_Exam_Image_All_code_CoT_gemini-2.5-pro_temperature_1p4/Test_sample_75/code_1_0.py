import math

def dms_to_decimal(degrees, minutes, seconds, direction):
    """Converts coordinates from DMS to decimal degrees."""
    decimal = degrees + minutes / 60 + seconds / 3600
    if direction in ['S', 'W']:
        decimal = -decimal
    return decimal

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculates the great-circle distance between two points on Earth.
    Returns the distance in miles.
    """
    R = 3958.8  # Earth radius in miles
    
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    distance = R * c
    return distance

# 1. Convert the rock carving's coordinates
carving_lat = dms_to_decimal(29, 3, 28.15, 'N')
carving_lon = dms_to_decimal(103, 48, 11.84, 'W')

# 2. Define coordinates for locations in the options
# Approx. center of the Chisos Mountains in Big Bend National Park
chisos_lat = 29.25
chisos_lon = -103.30
# Approx. coordinates for Lajitas, Texas
lajitas_lat = 29.263
lajitas_lon = -103.774

# 3. Calculate distances
dist_to_chisos = haversine_distance(carving_lat, carving_lon, chisos_lat, chisos_lon)
dist_to_lajitas = haversine_distance(carving_lat, carving_lon, lajitas_lat, lajitas_lon)

# 4. Print results and analysis
print("Step 1: Convert coordinates of the rock carving.")
print(f"Original coordinates: 29° 3' 28.15''N, 103° 48' 11.84''W")
print(f"Decimal coordinates: Lat {carving_lat:.4f}, Lon {carving_lon:.4f}\n")

print("Step 2: Evaluate geographical claims from the options.\n")

print("--- Analysis of Option B ---")
print("Claim: Carving depicts Chisos Mountains, about 20 miles north.")
print(f"Calculated distance to Chisos Mountains: {dist_to_chisos:.1f} miles.")
print("Conclusion: The claim of 'about 20 miles' is incorrect. The actual distance is over 30 miles.\n")

print("--- Analysis of Option D ---")
print("Claim: Carving depicts Bravo River near Lajitas, 10 miles northwest.")
print(f"Calculated distance to Lajitas: {dist_to_lajitas:.1f} miles.")
print("Conclusion: The claim of '10 miles' is incorrect. The actual distance is closer to 15 miles.\n")

print("Step 3: Final Conclusion.")
print("The geographical claims in options B and D are factually inaccurate based on distance calculations.")
print("Visually, the carving does not resemble a map but features symbolic shapes.")
print("The spiral carving strongly corresponds to the 'sipapu' symbol, a common motif in Southwestern indigenous rock art representing a portal or place of origin.")
print("Therefore, the most plausible interpretation is that the carving is symbolic, not geographical.")
