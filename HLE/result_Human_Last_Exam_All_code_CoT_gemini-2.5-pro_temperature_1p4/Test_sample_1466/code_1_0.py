import osmnx as ox
import networkx as nx
import warnings

# Suppress a common warning from osmnx for cleaner output
warnings.filterwarnings("ignore", category=UserWarning, message="The `inplace` parameter in pandas.Categorical.remove_unused_categories is deprecated")

def calculate_route_stats(graph, waypoints, route_name):
    """
    Calculates the total distance and estimated time for a route defined by waypoints.
    
    Args:
        graph (networkx.MultiDiGraph): The street network graph from osmnx.
        waypoints (list): A list of location names (strings) for the route.
        route_name (str): The name of the route (e.g., 'A').
        
    Returns:
        dict: A dictionary containing the route's name, distance, and time.
              Returns infinite distance/time on failure.
    """
    try:
        # Geocode the string waypoints to (latitude, longitude) coordinates
        coords = [ox.geocode(point) for point in waypoints]
        
        # Find the nearest nodes on the graph for each coordinate
        # We need to provide longitude first, then latitude to osmnx.nearest_nodes
        nodes = ox.nearest_nodes(graph, [c[1] for c in coords], [c[0] for c in coords])
        
        total_distance = 0
        # Calculate the shortest path length between each consecutive pair of nodes
        for i in range(len(nodes) - 1):
            start_node = nodes[i]
            end_node = nodes[i+1]
            # Use networkx's shortest_path_length with 'length' as the weight
            distance = nx.shortest_path_length(graph, start_node, end_node, weight='length')
            total_distance += distance
            
        # Average walking speed: 5 km/h = 5000 meters / 60 minutes = 83.33 m/min
        walking_speed_mpm = 83.33
        time_minutes = total_distance / walking_speed_mpm
        
        return {
            "name": route_name,
            "distance_m": total_distance,
            "time_min": time_minutes
        }
    except Exception as e:
        # Handle cases where a waypoint can't be found or a path doesn't exist
        print(f"Warning: Could not calculate route {route_name}. It might be invalid or disconnected. Reason: {e}")
        return {
            "name": route_name,
            "distance_m": float('inf'),
            "time_min": float('inf')
        }

def find_fastest_route():
    """Main function to define routes, calculate their stats, and find the fastest."""
    print("Downloading street network for the City of London...")
    # Define a bounding box around the area of interest for faster processing
    # Coords are (north, south, east, west)
    try:
        G = ox.graph_from_bbox(51.52, 51.51, -0.088, -0.105, network_type='walk')
    except Exception as e:
        print(f"Error: Could not download map data. Please check your internet connection. Details: {e}")
        return

    print("Analyzing routes...")

    # Define the key waypoints for each route to guide the path calculation
    routes_waypoints = {
        "A": [
            "Guildhall, London, UK",
            "Gresham Street and Foster Lane, London, UK",
            "Cheapside and New Change, London, UK",
            "St Paul's Cathedral, London, UK"
        ],
        "B": [
            "Guildhall, London, UK",
            "London Wall and Aldersgate Street, London, UK",
            "St Martin's Le Grand and Cheapside, London, UK",
            "St Paul's Cathedral, London, UK"
        ],
        "C": [
            "Guildhall, London, UK",
            "Lothbury and Princes Street, London, UK",
            "Queen Victoria Street and Cannon Street, London, UK",
            "St Paul's Cathedral, London, UK"
        ],
        "D": [
            "Guildhall, London, UK",
            "Cannon Street and St Swithin's Lane, London, UK",
            "Upper Thames Street, London, UK",
            "St Paul's Cathedral, London, UK"
        ],
        "E": [
            "Guildhall, London, UK",
            "Wood Street and London Wall, London, UK",
            "St Martin's Le Grand and Cheapside, London, UK",
            "St Paul's Cathedral, London, UK"
        ]
    }

    # Calculate statistics for all defined routes
    all_stats = []
    for name, waypoints in routes_waypoints.items():
        stats = calculate_route_stats(G, waypoints, name)
        all_stats.append(stats)
    
    # Filter out any routes that failed calculation
    valid_stats = [s for s in all_stats if s['time_min'] != float('inf')]
    
    if not valid_stats:
        print("Could not calculate statistics for any of the routes.")
        return

    # Find the route with the minimum walking time
    fastest_route = min(valid_stats, key=lambda x: x['time_min'])
    
    print("\n--- Route Analysis Results ---")
    print("Based on an average walking speed of 5 km/h (83.33 m/min):\n")
    
    walking_speed_mpm = 83.33
    # Sort by name for consistent output order
    for stats in sorted(valid_stats, key=lambda x: x['name']):
        dist = stats['distance_m']
        time = stats['time_min']
        # As requested, output each number in the final equation
        print(f"Route {stats['name']}: {dist:.0f} meters / {walking_speed_mpm:.2f} m/min = {time:.1f} minutes")
            
    print("\n-------------------------------------------")
    print(f"The fastest option is Route {fastest_route['name']}, taking approximately {fastest_route['time_min']:.1f} minutes.")
    print("-------------------------------------------")
    
    # Final answer in the specified format
    print(f"<<<{fastest_route['name']}>>>")

if __name__ == '__main__':
    find_fastest_route()