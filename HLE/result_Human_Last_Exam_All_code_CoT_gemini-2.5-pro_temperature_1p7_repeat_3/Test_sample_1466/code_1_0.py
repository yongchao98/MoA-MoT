import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance in kilometers between two points
    on the earth (specified in decimal degrees).
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    # Radius of earth in kilometers.
    r = 6371
    return c * r

# Define key waypoints for the routes in (latitude, longitude) format.
# These points define the general shape of each detour.
waypoints = {
    'start':        (51.5155, -0.0919),  # Guildhall
    'end':          (51.5138, -0.0984),  # St Paul's Cathedral
    'route_b_turn': (51.5176, -0.0952),  # Route B: Rotunda on London Wall (Northern detour)
    'route_c_turn': (51.5126, -0.0945),  # Route C: Queen Victoria St (Southern detour)
    'route_e_turn': (51.5180, -0.0963),  # Route E: Barbican area (Northern detour)
}

# A typical walking speed is about 5 km/h, which is 12 minutes per km.
WALKING_SPEED_MIN_PER_KM = 12

print("Estimating route times by calculating distances between key waypoints.\n")

# --- Route B Calculation ---
dist_b1 = haversine(waypoints['start'][0], waypoints['start'][1], waypoints['route_b_turn'][0], waypoints['route_b_turn'][1])
dist_b2 = haversine(waypoints['route_b_turn'][0], waypoints['route_b_turn'][1], waypoints['end'][0], waypoints['end'][1])
total_dist_b = dist_b1 + dist_b2
time_b = total_dist_b * WALKING_SPEED_MIN_PER_KM
print("--- Route B (Northern Detour via London Wall) ---")
print(f"Equation: ({dist_b1:.2f} km + {dist_b2:.2f} km) * {WALKING_SPEED_MIN_PER_KM} min/km = {time_b:.2f} minutes")
print(f"Total estimated distance: {total_dist_b:.2f} km\n")

# --- Route C Calculation ---
dist_c1 = haversine(waypoints['start'][0], waypoints['start'][1], waypoints['route_c_turn'][0], waypoints['route_c_turn'][1])
dist_c2 = haversine(waypoints['route_c_turn'][0], waypoints['route_c_turn'][1], waypoints['end'][0], waypoints['end'][1])
total_dist_c = dist_c1 + dist_c2
time_c = total_dist_c * WALKING_SPEED_MIN_PER_KM
print("--- Route C (Southern Detour via Queen Victoria St) ---")
print(f"Equation: ({dist_c1:.2f} km + {dist_c2:.2f} km) * {WALKING_SPEED_MIN_PER_KM} min/km = {time_c:.2f} minutes")
print(f"Total estimated distance: {total_dist_c:.2f} km\n")

# --- Route E Calculation ---
dist_e1 = haversine(waypoints['start'][0], waypoints['start'][1], waypoints['route_e_turn'][0], waypoints['route_e_turn'][1])
dist_e2 = haversine(waypoints['route_e_turn'][0], waypoints['route_e_turn'][1], waypoints['end'][0], waypoints['end'][1])
total_dist_e = dist_e1 + dist_e2
time_e = total_dist_e * WALKING_SPEED_MIN_PER_KM
print("--- Route E (Northern Detour via Barbican) ---")
print(f"Equation: ({dist_e1:.2f} km + {dist_e2:.2f} km) * {WALKING_SPEED_MIN_PER_KM} min/km = {time_e:.2f} minutes")
print(f"Total estimated distance: {total_dist_e:.2f} km\n")


# --- Conclusion ---
times = {
    "B": time_b,
    "C": time_c,
    "E": time_e,
}
fastest_route_letter = min(times, key=times.get)
print("="*30)
print(f"Based on the estimated distances, the fastest option is Route {fastest_route_letter}.")
print("="*30)
