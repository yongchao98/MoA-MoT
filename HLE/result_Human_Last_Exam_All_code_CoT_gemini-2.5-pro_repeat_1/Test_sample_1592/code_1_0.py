import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371
    return c * r

# Coordinates for Waskaganish, Quebec
waskaganish_lat = 51.49
waskaganish_lon = -78.76

# Coordinates for the closest points of neighboring provinces/territories
# Ontario: The point where the ON-QC border meets the James Bay coast
ontario_closest_lat = 51.47
ontario_closest_lon = -79.50

# Nunavut: The southeastern tip of Akimiski Island in James Bay
nunavut_closest_lat = 52.90
nunavut_closest_lon = -82.00


# Calculate distances
dist_to_ontario = haversine(waskaganish_lat, waskaganish_lon, ontario_closest_lat, ontario_closest_lon)
dist_to_nunavut = haversine(waskaganish_lat, waskaganish_lon, nunavut_closest_lat, nunavut_closest_lon)

# Compare distances and print the result
print("Finding the closest province or territory to Waskaganish (outside of Quebec)...")
print(f"Coordinates for Waskaganish, QC: ({waskaganish_lat}, {waskaganish_lon})")
print("-" * 30)

print(f"Distance to the closest point in Ontario ({ontario_closest_lat}, {ontario_closest_lon}): {dist_to_ontario:.2f} km")
print(f"Distance to the closest point in Nunavut ({nunavut_closest_lat}, {nunavut_closest_lon}): {dist_to_nunavut:.2f} km")
print("-" * 30)

if dist_to_ontario < dist_to_nunavut:
    closest_province = "Ontario"
    closest_dist = dist_to_ontario
else:
    closest_province = "Nunavut"
    closest_dist = dist_to_nunavut

print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest_province}.")
<<<Ontario>>>