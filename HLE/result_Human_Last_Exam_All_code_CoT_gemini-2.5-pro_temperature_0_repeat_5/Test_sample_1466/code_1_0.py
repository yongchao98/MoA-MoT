# First, please ensure you have the necessary libraries installed.
# You can install them by running the following command in your terminal:
# pip install osmnx numpy

import osmnx as ox
import networkx as nx
import numpy as np

def solve_routing_problem():
    """
    This function calculates the walking time for different routes in London
    to find the fastest option given a road closure.
    """
    print("Analyzing walking routes from Guildhall to St Paul's Cathedral...")
    print("Constraint: Cheapside road is closed between Grocers' Hall Court and Gutter Lane.")
    print("Route A is invalid as it uses the closed section of Cheapside.\n")

    # --- Step 1: Define a function to calculate route time ---
    def calculate_route_time(graph, waypoints_coords, avg_walk_speed_kmh=4.5):
        """
        Calculates the total length and time for a route defined by waypoints.
        """
        total_length_meters = 0
        try:
            # Find the path connecting all waypoints in sequence
            for i in range(len(waypoints_coords) - 1):
                start_point = waypoints_coords[i]
                end_point = waypoints_coords[i+1]
                
                # Find the nearest network nodes to our specified coordinates
                start_node = ox.nearest_nodes(graph, X=start_point[1], Y=start_point[0])
                end_node = ox.nearest_nodes(graph, X=end_point[1], Y=end_point[0])
                
                # Calculate the shortest path length between these two nodes
                path_length = nx.shortest_path_length(graph, source=start_node, target=end_node, weight='length')
                total_length_meters += path_length

            # Convert walking speed from km/h to meters/minute
            # Equation: speed_m_min = (speed_kmh * 1000) / 60
            avg_walk_speed_m_per_min = (avg_walk_speed_kmh * 1000) / 60
            
            # Calculate time in minutes
            # Equation: time_min = total_length_m / speed_m_min
            total_time_minutes = total_length_meters / avg_walk_speed_m_per_min
            
            return total_time_minutes, total_length_meters, avg_walk_speed_m_per_min
            
        except (nx.NetworkXNoPath, ValueError):
            # If a path is not found between waypoints
            return float('inf'), float('inf'), 0

    # --- Step 2: Get the street network and define waypoints ---
    print("Downloading London street network data... (This may take a moment)")
    # Define a central point and download the walking network graph
    graph_center_point = (51.5145, -0.095)
    G = ox.graph_from_point(graph_center_point, dist=1200, network_type='walk')
    print("Network data downloaded successfully.\n")

    # Define start, end, and route-specific waypoints
    guildhall_coords = (51.5155, -0.0920)
    st_pauls_coords = (51.5138, -0.0984)

    # Waypoints for Route B
    waypoints_b = [guildhall_coords, (51.5175, -0.0955), (51.5145, -0.0970), st_pauls_coords]
    # Waypoints for Route C
    waypoints_c = [guildhall_coords, (51.5140, -0.0890), (51.5125, -0.0940), st_pauls_coords]
    # Waypoints for Route D
    waypoints_d = [guildhall_coords, (51.5133, -0.0865), (51.5110, -0.0910), (51.5105, -0.0970), st_pauls_coords]
    # Waypoints for Route E
    waypoints_e = [guildhall_coords, (51.5175, -0.0940), (51.5175, -0.0955), st_pauls_coords]

    routes = {
        'B': waypoints_b,
        'C': waypoints_c,
        'D': waypoints_d,
        'E': waypoints_e
    }
    results = {}

    # --- Step 3: Calculate and display results for each route ---
    print("--- Calculating Route Times (at 4.5 km/h) ---\n")
    for name, waypoints in routes.items():
        time, length, speed = calculate_route_time(G, waypoints)
        if np.isinf(time):
            print(f"Route {name}: Could not find a valid path.")
        else:
            # Outputting the numbers used in the final calculation for this route
            print(f"Route {name}:")
            print(f"  - Total Distance: {length:.0f} meters")
            print(f"  - Walking Speed: {speed:.0f} meters/minute")
            print(f"  - Calculation: {length:.0f} (meters) / {speed:.0f} (m/min) = {time:.2f} minutes")
            print(f"  - Estimated Time: {time:.2f} minutes\n")
        results[name] = time

    # --- Step 4: Determine the fastest route ---
    # Filter out any routes that might have failed
    valid_results = {k: v for k, v in results.items() if not np.isinf(v)}
    if not valid_results:
        print("Could not determine the fastest route.")
        return

    fastest_route = min(valid_results, key=valid_results.get)
    fastest_time = valid_results[fastest_route]

    print("--- Conclusion ---")
    print(f"Comparing the estimated times, the fastest option is Route {fastest_route}.")
    
    # --- Final Answer ---
    print(f"\n<<<{fastest_route}>>>")

# Run the main function
solve_routing_problem()