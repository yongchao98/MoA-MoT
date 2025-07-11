import geopy.distance

# Coordinates for the locations
# Waskaganish, Quebec
waskaganish_lat, waskaganish_lon = 51.487, -78.758

# The Ontario-Quebec border meets James Bay at approximately longitude 79°31′W.
# The closest point on the Ontario border is on the James Bay coast.
ontario_border_lat, ontario_border_lon = 51.467, -79.517

# The closest part of the territory of Nunavut is Akimiski Island in James Bay.
# We'll use a point on the southeastern coast of the island.
nunavut_island_lat, nunavut_island_lon = 52.933, -79.450

# Define the coordinate tuples
waskaganish_coords = (waskaganish_lat, waskaganish_lon)
ontario_coords = (ontario_border_lat, ontario_border_lon)
nunavut_coords = (nunavut_island_lat, nunavut_island_lon)

# Calculate the geodesic distance using the geopy library
distance_to_ontario_km = geopy.distance.geodesic(waskaganish_coords, ontario_coords).km
distance_to_nunavut_km = geopy.distance.geodesic(waskaganish_coords, nunavut_coords).km

# Print the information used for the calculation
print(f"Finding the closest province or territory to Waskaganish, QC ({waskaganish_lat}, {waskaganish_lon})")
print("-" * 30)

# Print the calculation for Ontario
print("Candidate: Ontario")
print(f"Closest Point (approx.): ({ontario_border_lat}, {ontario_border_lon})")
print(f"Distance from Waskaganish to Ontario border: {distance_to_ontario_km:.2f} km")
print("-" * 30)

# Print the calculation for Nunavut
print("Candidate: Nunavut")
print(f"Closest Point (Akimiski Island, approx.): ({nunavut_island_lat}, {nunavut_island_lon})")
print(f"Distance from Waskaganish to Nunavut: {distance_to_nunavut_km:.2f} km")
print("-" * 30)

# Determine and print the final answer
if distance_to_ontario_km < distance_to_nunavut_km:
    closest_place = "Ontario"
    closest_distance = distance_to_ontario_km
else:
    closest_place = "Nunavut"
    closest_distance = distance_to_nunavut_km

print(f"Conclusion: The closest province or territory to Waskaganish is {closest_place}, at an approximate distance of {closest_distance:.2f} km.")

<<<Ontario>>>