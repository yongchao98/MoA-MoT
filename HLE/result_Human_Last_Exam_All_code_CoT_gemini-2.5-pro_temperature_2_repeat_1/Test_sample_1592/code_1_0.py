import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    r = 6371  # Radius of earth in kilometers.
    return c * r

# Coordinates of the locations
waskaganish_qc = {"name": "Waskaganish, QC", "lat": 51.48, "lon": -78.76}
ontario_border = {"name": "Ontario Border", "lat": 51.48, "lon": -79.52} # Approx. location of the ON-QC border at the same latitude
nunavut_island = {"name": "Akimiski Island, NU", "lat": 52.99, "lon": -81.33} # Southernmost tip of the closest island in Nunavut

# Calculate distances
dist_to_ontario = haversine_distance(waskaganish_qc["lat"], waskaganish_qc["lon"], ontario_border["lat"], ontario_border["lon"])
dist_to_nunavut = haversine_distance(waskaganish_qc["lat"], waskaganish_qc["lon"], nunavut_island["lat"], nunavut_island["lon"])

print(f"Calculating distances from {waskaganish_qc['name']}...")
print(f"Distance to the nearest point in Ontario: {dist_to_ontario:.2f} km")
print(f"Distance to the nearest point in Nunavut (Akimiski Island): {dist_to_nunavut:.2f} km")
print("-" * 30)

# Determine the closest province/territory
if dist_to_ontario < dist_to_nunavut:
    closest_place = "Ontario"
else:
    closest_place = "Nunavut"

print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest_place}.")

<<<Ontario>>>