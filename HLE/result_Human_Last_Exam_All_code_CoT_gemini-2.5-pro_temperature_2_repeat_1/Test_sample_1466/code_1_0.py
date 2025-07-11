import geopy.distance

# Define the waypoints for Route E based on the directions provided.
# Each waypoint is a tuple: (Description, Latitude, Longitude)
route_e_waypoints = [
    ("Start: Guildhall", 51.5155, -0.0919),
    ("Turn towards Wood St", 51.5160, -0.0933),
    ("Barbican Highwalk / Gilbert Bridge", 51.5175, -0.0940),
    ("Junction with Aldersgate St", 51.5178, -0.0960),
    ("St Martin's Le Grand", 51.5155, -0.0967),
    ("Junction near Cheapside", 51.5143, -0.0975),
    ("End: St Paul's Cathedral", 51.5138, -0.0984)
]

total_distance_meters = 0
distance_segments = []

# Calculate the distance for each segment of the journey
for i in range(len(route_e_waypoints) - 1):
    point1 = (route_e_waypoints[i][1], route_e_waypoints[i][2])
    point2 = (route_e_waypoints[i+1][1], route_e_waypoints[i+1][2])
    segment_distance = geopy.distance.great_circle(point1, point2).meters
    distance_segments.append(segment_distance)
    total_distance_meters += segment_distance

# Estimate walking time assuming an average speed of 5 km/h (83.33 meters/minute)
walking_speed_m_per_min = 83.33
walking_time_minutes = total_distance_meters / walking_speed_m_per_min

# Output the calculations as requested
print("Calculation for estimated walking time of Route E:")
equation_parts = [f"{dist:.0f}m" for dist in distance_segments]
print(" + ".join(equation_parts) + f" = {total_distance_meters:.0f} meters")
print(f"\n{total_distance_meters:.0f} meters / {walking_speed_m_per_min:.2f} meters_per_minute = {walking_time_minutes:.1f} minutes")
print(f"\nThe estimated fastest walking time is approximately {round(walking_time_minutes)} minutes via Route E.")
