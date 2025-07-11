import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance in kilometers between two points
    on the earth (specified in decimal degrees).
    """
    # Earth radius in kilometers
    R = 6371.0

    # Convert decimal degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Haversine formula
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = R * c
    return distance

# Coordinates for the locations
# 1. Waskaganish, Quebec
waskaganish_lat, waskaganish_lon = 51.1944, -78.7592

# 2. Closest point in Ontario (ON/QC border on James Bay coast)
ontario_lat, ontario_lon = 51.47, -79.50

# 3. Closest point in Nunavut (SE tip of Akimiski Island)
nunavut_lat, nunavut_lon = 52.90, -81.30

# Calculate distances
dist_to_ontario = haversine_distance(waskaganish_lat, waskaganish_lon, ontario_lat, ontario_lon)
dist_to_nunavut = haversine_distance(waskaganish_lat, waskaganish_lon, nunavut_lat, nunavut_lon)

# Print the comparison
print(f"Finding the closest province or territory to Waskaganish, QC ({waskaganish_lat:.2f} N, {abs(waskaganish_lon):.2f} W)...")
print("-" * 20)

print("Candidate 1: Ontario")
print(f"Closest point in Ontario is at {ontario_lat:.2f} N, {abs(ontario_lon):.2f} W.")
print(f"Equation: distance(Waskaganish({waskaganish_lat}, {waskaganish_lon}), Ontario({ontario_lat}, {ontario_lon}))")
print(f"Result: The distance is approximately {dist_to_ontario:.1f} km.")
print("-" * 20)

print("Candidate 2: Nunavut")
print(f"Closest point in Nunavut is at {nunavut_lat:.2f} N, {abs(nunavut_lon):.2f} W.")
print(f"Equation: distance(Waskaganish({waskaganish_lat}, {waskaganish_lon}), Nunavut({nunavut_lat}, {nunavut_lon}))")
print(f"Result: The distance is approximately {dist_to_nunavut:.1f} km.")
print("-" * 20)

# Determine the closest one
if dist_to_ontario < dist_to_nunavut:
    closest_location = "Ontario"
else:
    closest_location = "Nunavut"

print(f"Comparing the distances, {closest_location} is closer to Waskaganish than the other candidate.")

# Final answer in the required format
# print(f"<<<{closest_location}>>>")