# First, you may need to install the geopy library:
# pip install geopy

from geopy.distance import great_circle

# Coordinates for the Waskaganish gathering place in Quebec
waskaganish_coords = (51.1939, -78.7622)

# 1. ONTARIO
# The Quebec-Ontario border south of James Bay is along the meridian 79° 31' 0" W.
# 31 minutes = 31/60 = 0.5167 degrees. So the longitude is -79.5167.
# We'll find the distance to a point on that border at the same latitude as Waskaganish.
ontario_border_coords = (waskaganish_coords[0], -79.5167)
dist_to_ontario = great_circle(waskaganish_coords, ontario_border_coords).km

# 2. NUNAVUT
# The closest major landmass in Nunavut is Akimiski Island in James Bay.
# We'll use the coordinates of a point on the southeastern coast of the island.
nunavut_akimiksi_island_coords = (52.25, -79.5)
dist_to_nunavut = great_circle(waskaganish_coords, nunavut_akimiksi_island_coords).km

# Compare the distances and print the result
print(f"The Waskaganish gathering place is located at {waskaganish_coords[0]:.2f}° N, {abs(waskaganish_coords[1]):.2f}° W.")
print("-" * 30)
print(f"Distance to the Ontario border: {dist_to_ontario:.2f} km")
print(f"Distance to the Nunavut border (Akimiski Island): {dist_to_nunavut:.2f} km")
print("-" * 30)

if dist_to_ontario < dist_to_nunavut:
    closest = "Ontario"
    closest_dist = dist_to_ontario
else:
    closest = "Nunavut"
    closest_dist = dist_to_nunavut

print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest}.")
<<<Ontario>>>