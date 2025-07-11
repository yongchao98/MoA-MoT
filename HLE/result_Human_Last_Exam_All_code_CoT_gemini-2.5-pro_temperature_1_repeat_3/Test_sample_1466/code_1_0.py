import osmnx as ox
import networkx as nx
import warnings

# Suppress a common warning from osmnx about query strings
warnings.filterwarnings("ignore", category=UserWarning, module='osmnx')

def calculate_route_stats(graph, waypoints_coords, walking_speed_mps=1.4):
    """
    Calculates the total distance and time for a route defined by a list of coordinates.
    A route is a sequence of waypoints, and the total distance is the sum of the
    shortest paths between consecutive waypoints.
    """
    total_distance = 0
    
    try:
        # Get the nearest graph nodes for all waypoint coordinates
        nodes = ox.nearest_nodes(graph, [coord[1] for coord in waypoints_coords], [coord[0] for coord in waypoints_coords])
    except Exception:
        # Return infinite distance/time if nodes can't be found
        return float('inf'), float('inf')

    # Calculate the path length between each consecutive pair of nodes
    for i in range(len(nodes) - 1):
        try:
            path_length = nx.shortest_path_length(graph, source=nodes[i], target=nodes[i+1], weight='length')
            total_distance += path_length
        except nx.NetworkXNoPath:
            # If no path exists between two waypoints, the route is invalid
            return float('inf'), float('inf')
            
    # Convert distance (meters) to time (minutes)
    total_time_seconds = total_distance / walking_speed_mps
    total_time_minutes = total_time_seconds / 60
    
    return total_distance, total_time_minutes

def solve_fastest_route():
    """
    Main function to find the fastest walking route among the given options.
    """
    # This script requires the osmnx and networkx libraries.
    # You can install them with: pip install osmnx networkx
    print("Initializing...")

    # Define the start, end, and key intermediate points for each route
    # These waypoints are chosen to accurately represent each route option.
    start_point = "Guildhall, London, UK"
    end_point = "St Paul's Cathedral, London, UK"
    
    routes_definitions = {
        'A': [start_point, "Foster Lane, London, UK", end_point],
        'B': [start_point, "Rotunda, Aldersgate St, London, UK", end_point],
        'C': [start_point, "Mansion House Station, London, UK", end_point],
        'D': [start_point, "Cannon Street Station, London, UK", "High Timber Street, London, UK", end_point],
        'E': [start_point, "Barbican Centre, London, UK", end_point]
    }

    # Geocode all unique locations to get (latitude, longitude) coordinates
    # This is done in a batch to be efficient.
    unique_locations = set()
    for path in routes_definitions.values():
        unique_locations.update(path)
    
    print("Fetching coordinates for key locations...")
    try:
        # The 'pause' parameter respects the Nominatim API usage policy.
        coords_list = ox.geocode_many(list(unique_locations), pause=1)
        coords_dict = dict(zip(list(unique_locations), coords_list))
    except Exception as e:
        print(f"Error: Could not fetch coordinates. Please check your internet connection. Details: {e}")
        return

    # Download the walking network graph for the area from OpenStreetMap
    print("Downloading map data for the area (this may take a moment)...")
    try:
        center_point = coords_dict[start_point]
        graph = ox.graph_from_point(center_point, dist=1500, network_type='walk')
    except Exception as e:
        print(f"Error: Could not download map data. Please check your internet connection. Details: {e}")
        return

    # Calculate and store the results for each route
    print("\nCalculating walking time for each route...")
    results = {}
    for route_id, path in routes_definitions.items():
        path_coords = [coords_dict[loc] for loc in path if loc in coords_dict and coords_dict[loc] is not None]
        if len(path_coords) == len(path):
            distance, time = calculate_route_stats(graph, path_coords)
            results[route_id] = (distance, time)
        else:
            results[route_id] = (float('inf'), float('inf'))

    print("\n--- Route Analysis Complete ---\n")

    # Print the final equations (results) for each route
    for route_id, (distance, time) in sorted(results.items()):
        if time == float('inf'):
            print(f"Route {route_id}: Could not be calculated (invalid path).")
        else:
            # Each number in the final equation: Route, Distance, Time
            print(f"Route {route_id}: Total Distance = {distance:.0f} meters; Estimated Time = {time:.1f} minutes")

    # Find and announce the fastest route
    if all(time == float('inf') for _, time in results.values()):
        print("\nCould not determine the fastest route.")
    else:
        fastest_route_id = min(results, key=lambda r: results[r][1])
        min_time = results[fastest_route_id][1]
        print("\n-------------------------------------------------")
        print(f"The fastest route is Route {fastest_route_id} with an estimated time of {min_time:.1f} minutes.")
        print("-------------------------------------------------")

if __name__ == "__main__":
    solve_fastest_route()