# First, ensure you have the necessary libraries installed:
# pip install osmnx networkx

import osmnx as ox
import networkx as nx

def calculate_route_metrics(graph, waypoints, speed_kph=5.0):
    """
    Calculates the total distance and walking time for a route defined by a series of waypoints.

    Args:
        graph (networkx.MultiDiGraph): The street network graph from osmnx.
        waypoints (list): A list of (latitude, longitude) tuples representing the route.
        speed_kph (float): Average walking speed in kilometers per hour.

    Returns:
        tuple: A tuple containing the total distance in meters and the total time in minutes.
               Returns (float('inf'), float('inf')) if a path is not found.
    """
    total_length_meters = 0
    # Convert speed from km/h to m/min for time calculation
    speed_m_per_min = speed_kph * 1000 / 60

    try:
        # Iterate through pairs of consecutive waypoints to build the full route
        for i in range(len(waypoints) - 1):
            start_point = waypoints[i]
            end_point = waypoints[i+1]

            # Find the nearest graph nodes to the geographic coordinates
            start_node = ox.distance.nearest_nodes(graph, X=start_point[1], Y=start_point[0])
            end_node = ox.distance.nearest_nodes(graph, X=end_point[1], Y=end_point[0])

            # Calculate the shortest path length for this leg of the journey
            leg_length = nx.shortest_path_length(graph, source=start_node, target=end_node, weight='length')
            total_length_meters += leg_length

        # Calculate total walking time in minutes
        total_time_minutes = total_length_meters / speed_m_per_min
        return total_length_meters, total_time_minutes
    except nx.NetworkXNoPath:
        # This can happen if waypoints are in disconnected parts of the graph
        return float('inf'), float('inf')

# --- Main script execution ---

# Define the start (Guildhall) and end (St Paul's Cathedral) coordinates
guildhall_coord = (51.5156, -0.0920)
st_pauls_coord = (51.5138, -0.0984)

# A central point for downloading the map area
center_point = (51.514, -0.094)
print("Downloading London street network data for walking paths...")
# Download the walkable street network graph
G = ox.graph_from_point(center_point, dist=1500, network_type='walk')
print("Download complete.\n")

# Define waypoints for each valid route (B, C, D, E) based on their descriptions.
# Route A is excluded as it is invalid.
routes_info = {
    'B': {
        "name": "Northern Route via London Wall",
        "waypoints": [
            guildhall_coord,
            (51.5173, -0.0958),  # Aldersgate St/London Wall Rotunda
            st_pauls_coord
        ]
    },
    'C': {
        "name": "Eastern/Southern Route via Queen Victoria St",
        "waypoints": [
            guildhall_coord,
            (51.5120, -0.0933),  # Queen Victoria St / Cannon St junction
            st_pauls_coord
        ]
    },
    'D': {
        "name": "Far Southern Route via Upper Thames St",
        "waypoints": [
            guildhall_coord,
            (51.5121, -0.0886),  # Cannon St / St Swithin's Ln
            (51.5103, -0.0967),  # High Timber St area
            st_pauls_coord
        ]
    },
    'E': {
        "name": "Far Northern Route via Barbican",
        "waypoints": [
            guildhall_coord,
            (51.5179, -0.0945),  # St Giles Terrace near Barbican
            st_pauls_coord
        ]
    }
}

# Calculate and store the results for each route
results = {}
print("Calculating times for each possible route...\n")

for key, data in routes_info.items():
    distance, time = calculate_route_metrics(G, data['waypoints'])
    results[key] = {'time': time, 'distance': distance, 'name': data['name']}
    print(f"--- Route {key}: {data['name']} ---")
    if time == float('inf'):
        print("  Could not find a valid path for this route.")
    else:
        # Output the numbers used in the calculation
        print(f"  Total Distance: {distance:.0f} meters")
        print(f"  Estimated Walking Time: {time:.1f} minutes")
    print()

# Determine the fastest route from the calculated results
if not any(results[r]['time'] != float('inf') for r in results):
    print("Could not determine the fastest route.")
else:
    fastest_route_key = min(results, key=lambda k: results[k]['time'])
    fastest_route_data = results[fastest_route_key]

    print("=============================================")
    print("                CONCLUSION")
    print("=============================================")
    print(f"The fastest walking route is Route {fastest_route_key}.")
    print(f"It covers a distance of {fastest_route_data['distance']:.0f} meters and takes approximately {fastest_route_data['time']:.1f} minutes.")
