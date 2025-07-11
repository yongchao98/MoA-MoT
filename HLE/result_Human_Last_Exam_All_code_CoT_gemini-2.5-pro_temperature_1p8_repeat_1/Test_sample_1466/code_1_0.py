import osmnx as ox
import networkx as nx

def calculate_route_time(graph, waypoints, walking_speed_kph=5.0):
    """Calculates the total distance and time for a route defined by waypoints."""
    total_distance = 0
    try:
        # Find the nearest graph nodes for the waypoint coordinates
        nodes = ox.distance.nearest_nodes(graph, [p[1] for p in waypoints], [p[0] for p in waypoints])
    except Exception:
        # If a point is outside the graph, the route is invalid for this map
        return float('inf'), float('inf')

    # Calculate the path length between consecutive waypoints
    for i in range(len(nodes) - 1):
        try:
            length = nx.shortest_path_length(graph, nodes[i], nodes[i+1], weight='length')
            total_distance += length
        except nx.NetworkXNoPath:
            # If no path exists between two waypoints, the route is invalid.
            return float('inf'), float('inf')

    # Calculate walking time in minutes from distance in meters
    time_minutes = (total_distance / 1000) / walking_speed_kph * 60
    return total_distance, time_minutes

def find_fastest_walk():
    """
    This function calculates the fastest walking route from a set of options
    by modeling the routes on a real city map.
    """
    # Define Start and End points
    guildhall = (51.5155, -0.0920)
    st_pauls = (51.5138, -0.0984)

    # Define a bounding box around the area to download map data
    north, south, east, west = 51.519, 51.512, -0.088, -0.099

    try:
        print("Downloading map data for Central London... (This may take a moment)")
        # Download the walking network graph for the specified area
        G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
        print("Map data downloaded successfully.\n")
    except Exception as e:
        print(f"Could not download map data. Please check your internet connection.")
        print(f"Error: {e}")
        return

    print("Calculating and comparing route times...\n")
    
    # --- Define Route Waypoints ---
    # Waypoints are used to guide the pathfinding along the described routes.

    # Route A (Short detour via Gresham St and Foster Ln)
    waypoints_a = [guildhall, (51.5150, -0.0960), st_pauls]
    dist_a, time_a = calculate_route_time(G, waypoints_a)

    # Route B (Long northern detour via London Wall)
    waypoints_b = [guildhall, (51.5173, -0.0917), (51.5181, -0.0969), st_pauls]
    dist_b, time_b = calculate_route_time(G, waypoints_b)

    # Route C (Southern detour via Queen Victoria St)
    waypoints_c = [guildhall, (51.5144, -0.0895), (51.5131, -0.0906), st_pauls]
    dist_c, time_c = calculate_route_time(G, waypoints_c)
    
    # --- Print Results with Calculation Steps ---
    print("--- Analysis of Potential Routes (assuming 5 km/h walking speed) ---\n")
    
    # Route A Results
    print("--- Route A (Short Detour via Foster Ln) ---")
    print(f"Distance: {dist_a:.0f} meters")
    print(f"Calculation: ({dist_a:.0f} / 1000) / 5 * 60 = {time_a:.1f} minutes")
    print(f"Estimated Time: {time_a:.1f} minutes\n")

    # Route B Results
    print("--- Route B (Northern Detour) ---")
    print(f"Distance: {dist_b:.0f} meters")
    print(f"Calculation: ({dist_b:.0f} / 1000) / 5 * 60 = {time_b:.1f} minutes")
    print(f"Estimated Time: {time_b:.1f} minutes\n")

    # Route C Results
    print("--- Route C (Southern Detour) ---")
    print(f"Distance: {dist_c:.0f} meters")
    print(f"Calculation: ({dist_c:.0f} / 1000) / 5 * 60 = {time_c:.1f} minutes")
    print(f"Estimated Time: {time_c:.1f} minutes\n")

    # --- Conclusion ---
    print("--- Comparison ---")
    times = {'A': time_a, 'B': time_b, 'C': time_c}
    valid_times = {k: v for k, v in times.items() if v != float('inf')}
    
    if not valid_times:
        print("Could not find a valid path for any route.")
    else:
        fastest_route_letter = min(valid_times, key=valid_times.get)
        print(f"Based on the calculations, Route {fastest_route_letter} is the fastest option.")


if __name__ == '__main__':
    # To run this script, you need to have osmnx and its dependencies installed:
    # pip install osmnx
    find_fastest_walk()