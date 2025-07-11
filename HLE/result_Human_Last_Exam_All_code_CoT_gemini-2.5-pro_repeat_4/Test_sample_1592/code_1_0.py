import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
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

# Coordinates of Waskaganish, Quebec
waskaganish_lat, waskaganish_lon = 51.483, -78.767

# Coordinates of the closest point in Ontario (provincial border on James Bay coast)
ontario_point_lat, ontario_point_lon = 51.483, -79.517

# Coordinates of the closest point in Nunavut (southeastern coast of Akimiski Island)
nunavut_point_lat, nunavut_point_lon = 52.8, -80.8

# Calculate distances
dist_to_ontario = haversine_distance(waskaganish_lat, waskaganish_lon, ontario_point_lat, ontario_point_lon)
dist_to_nunavut = haversine_distance(waskaganish_lat, waskaganish_lon, nunavut_point_lat, nunavut_point_lon)

print(f"The Waskaganish gathering place is located at latitude {waskaganish_lat}, longitude {waskaganish_lon}.")
print(f"The closest point in Ontario is at latitude {ontario_point_lat}, longitude {ontario_point_lon}.")
print(f"The closest point in Nunavut (Akimiski Island) is at latitude {nunavut_point_lat}, longitude {nunavut_point_lon}.")
print("-" * 30)
print(f"Distance to Ontario border: {dist_to_ontario:.2f} km")
print(f"Distance to Nunavut (Akimiski Island): {dist_to_nunavut:.2f} km")
print("-" * 30)

# Determine the closest province/territory
if dist_to_ontario < dist_to_nunavut:
    closest_location = "Ontario"
else:
    closest_location = "Nunavut"

print(f"The closest province or territory to the Waskaganish gathering place, outside of Quebec, is {closest_location}.")

<<<Ontario>>>