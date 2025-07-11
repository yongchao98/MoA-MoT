from geopy.distance import geodesic
import operator

# Geographic coordinates for the locations
# Waskaganish, QC is our starting point
waskaganish_coords = (51.1947, -78.7611)

# Representative points for neighboring provinces/territories
locations = {
    "Ontario": (51.2731, -80.6471),        # Moosonee, ON
    "Nunavut": (52.41, -80.85),            # Southernmost tip of Akimiski Island, NU
    "Newfoundland and Labrador": (52.9366, -66.8788) # Wabush, NL
}

# Dictionary to store the calculated distances
distances = {}

print(f"Calculating distances from Waskaganish, QC {waskaganish_coords} to neighboring regions:\n")

# Calculate the distance to each location
for province, coords in locations.items():
    distance = geodesic(waskaganish_coords, coords).km
    distances[province] = distance
    print(f"Distance to {province} (via {coords}): {distance:.2f} km")

# Find the location with the minimum distance
closest_province = min(distances.items(), key=operator.itemgetter(1))

print(f"\nThe closest province or territory to Waskaganish, outside of Quebec, is {closest_province[0]}.")

<<<Ontario>>>