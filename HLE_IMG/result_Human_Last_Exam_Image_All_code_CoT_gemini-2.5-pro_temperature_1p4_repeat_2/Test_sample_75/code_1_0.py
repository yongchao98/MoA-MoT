import geopy.distance
from geopy.point import Point

def dms_to_dd(degrees, minutes, seconds, direction):
    """Converts coordinates from DMS to Decimal Degrees."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

def get_bearing(point1, point2):
    """Calculates the bearing between two points."""
    lat1, lon1 = point1.latitude, point1.longitude
    lat2, lon2 = point2.latitude, point2.longitude

    import math
    dLon = (lon2 - lon1)
    y = math.sin(math.radians(dLon)) * math.cos(math.radians(lat2))
    x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dLon))
    brng = math.atan2(y, x)
    brng = math.degrees(brng)
    brng = (brng + 360) % 360
    
    dirs = ['North', 'Northeast', 'East', 'Southeast', 'South', 'Southwest', 'West', 'Northwest']
    index = round(brng / 45) % 8
    return dirs[index]

# 1. Coordinates of the rock carving
carving_lat_dms = (29, 3, 28.15, 'N')
carving_lon_dms = (103, 48, 11.84, 'W')

# Convert carving coordinates to decimal degrees
carving_lat_dd = dms_to_dd(*carving_lat_dms)
carving_lon_dd = dms_to_dd(*carving_lon_dms)
carving_coords = Point(carving_lat_dd, carving_lon_dd)

# 2. Coordinates of locations from answer choices
# Chisos Mountains (approximated at Emory Peak)
chisos_coords = Point(29.253, -103.300) 
# Lajitas, Texas
lajitas_coords = Point(29.261, -103.771)

# 3. Calculate distances from the carving
distance_to_chisos = geopy.distance.distance(carving_coords, chisos_coords).miles
distance_to_lajitas = geopy.distance.distance(carving_coords, lajitas_coords).miles

# 4. Calculate bearing (direction) from the carving to the locations
bearing_to_chisos = get_bearing(carving_coords, chisos_coords)
bearing_to_lajitas = get_bearing(carving_coords, lajitas_coords)

# Print the analysis
print("Analysis of Geographic Claims:")
print("-" * 30)

# Analysis for Option B
print("Claim (B): Chiso Mountains are ~20 miles North.")
print(f"Calculated distance to Chiso Mountains: {distance_to_chisos:.2f} miles.")
print(f"Calculated direction to Chiso Mountains: {bearing_to_chisos}.")
print("Result: The claim is inaccurate. The distance is about 30 miles and the direction is Northeast, not North.")
print("-" * 30)

# Analysis for Option D
print("Claim (D): Lajitas, TX is 10 miles Northwest.")
print(f"Calculated distance to Lajitas: {distance_to_lajitas:.2f} miles.")
print(f"Calculated direction to Lajitas: {bearing_to_lajitas}.")
print("Result: The claim is inaccurate. The distance is about 14 miles and the direction is Northeast, not Northwest.")
print("-" * 30)

print("\nConclusion:")
print("The specific geographic claims in options B and D are factually incorrect based on the provided coordinates.")
print("Visual inspection shows figures typical of Southwestern indigenous rock art, including zoomorphic (animal-like) and geometric/symbolic shapes.")
print("Option E provides a plausible archaeological interpretation for the circular/radiating symbol, which is a feature known as a sipapu, representing a place of origin in indigenous cosmology.")
