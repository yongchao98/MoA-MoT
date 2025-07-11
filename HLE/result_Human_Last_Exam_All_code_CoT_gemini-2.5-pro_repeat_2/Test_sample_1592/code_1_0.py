import math

def calculate_haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculates the distance between two points on Earth using the Haversine formula.
    """
    R = 6371  # Radius of Earth in kilometers
    
    # Convert latitude and longitude from degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    
    # Difference in coordinates
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    
    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    distance = R * c
    return distance

# Coordinates for the locations
# 1. Waskaganish Gathering Place, Quebec
waskaganish_lat, waskaganish_lon = 51.1944, -78.7606

# 2. Closest point in Ontario (on the James Bay coast, just across the provincial border)
ontario_border_lat, ontario_border_lon = 51.20, -79.52

# 3. Closest point in Nunavut (South Twin Island, the southernmost island of Nunavut in James Bay)
nunavut_island_lat, nunavut_island_lon = 53.13, -79.88

# Calculate distances
distance_to_ontario = calculate_haversine_distance(waskaganish_lat, waskaganish_lon, ontario_border_lat, ontario_border_lon)
distance_to_nunavut = calculate_haversine_distance(waskaganish_lat, waskaganish_lon, nunavut_island_lat, nunavut_island_lon)

# Print the results
print(f"Finding the closest province or territory to Waskaganish, QC (outside of Quebec)...")
print("-" * 50)
print(f"Coordinates for Waskaganish, QC: Latitude={waskaganish_lat}, Longitude={waskaganish_lon}")
print(f"Coordinates for closest point in Ontario: Latitude={ontario_border_lat}, Longitude={ontario_border_lon}")
print(f"Coordinates for closest point in Nunavut: Latitude={nunavut_island_lat}, Longitude={nunavut_island_lon}")
print("-" * 50)
print(f"Distance from Waskaganish to the Ontario border: {distance_to_ontario:.2f} km")
print(f"Distance from Waskaganish to Nunavut (South Twin Island): {distance_to_nunavut:.2f} km")
print("-" * 50)

# Determine and print the closest one
if distance_to_ontario < distance_to_nunavut:
    closest_place = "Ontario"
else:
    closest_place = "Nunavut"

print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest_place}.")

<<<Ontario>>>