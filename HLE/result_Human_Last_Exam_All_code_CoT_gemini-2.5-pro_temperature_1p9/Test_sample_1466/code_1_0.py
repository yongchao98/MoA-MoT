import osmnx as ox
import networkx as nx
import sys

# Suppress osmnx informational messages for a cleaner output.
# osmnx caches map data locally to speed up subsequent runs.
ox.settings.log_console = False
ox.settings.use_cache = True

def solve_walking_route():
    """
    This function analyzes several walking routes from Guildhall to St Paul's Cathedral,
    calculating the distance and estimated time for each to find the fastest option
    given a road closure on Cheapside.
    """
    print("Analyzing walking routes from Guildhall to St Paul's Cathedral...")
    print("A section of Cheapside is closed, so we must find the best detour.")
    print("The fastest route will be the one with the shortest total distance.\n")

    # Define a central point and radius for the map download to ensure all routes are covered.
    map_center_point = (51.515, -0.095)
    try:
        # Download the walkable street network from OpenStreetMap
        G = ox.graph_from_point(map_center_point, dist=1500, network_type='walk')
    except Exception as e:
        print(f"Failed to download map data. Please check your internet connection. Error: {e}", file=sys.stderr)
        return

    # Define the waypoints for each route based on the provided answer choices.
    # These are representative junctions or landmarks to trace each path.
    routes_waypoints = {
        'A': [
            "Guildhall, London, UK",
            "intersection of Gresham Street and Foster Lane, London, UK",
            "intersection of Foster Lane and Cheapside, London, UK",
            "St Paul's Cathedral, London, UK"
        ],
        'B': [
            "Guildhall, London, UK",
            "Museum of London Roundabout, London, UK", # Represents the 'Rotunda' at London Wall
            "intersection of St. Martin's Le Grand and Cheapside, London, UK",
            "St Paul's Cathedral, London, UK"
        ],
        'C': [
            "Guildhall, London, UK",
            "intersection of Gresham Street and Princes Street, London, UK",
            "Mansion House, London, UK", # A landmark on Queen Victoria Street
            "St Paul's Cathedral, London, UK"
        ],
        'D': [
            "Guildhall, London, UK",
            "intersection of Lothbury and Princes Street, London, UK",
            "intersection of Cannon Street and St Swithin's Lane, London, UK",
            "Millennium Bridge, London, UK", # Represents the path near the river
            "St Paul's Cathedral, London, UK"
        ],
        'E': [
            "Guildhall, London, UK",
            "Barbican Centre, London, UK", # A landmark to the north
            "intersection of Aldersgate Street and London Wall, London, UK",
            "St Paul's Cathedral, London, UK"
        ]
    }

    # Average walking speed: 4.5 km/h, which is 75 meters per minute.
    walking_speed_mpm = 75

    results = {}

    print("Calculating distances and times for each route...\n")

    # Geocode waypoints and calculate path lengths
    for route_id, waypoints in routes_waypoints.items():
        try:
            # Get the graph nodes closest to each waypoint's geographic coordinates
            waypoint_nodes = [ox.nearest_nodes(G, ox.geocode(wp)[1], ox.geocode(wp)[0]) for wp in waypoints]
            
            total_length = 0
            # Calculate the shortest path length between each consecutive waypoint
            for i in range(len(waypoint_nodes) - 1):
                source_node = waypoint_nodes[i]
                target_node = waypoint_nodes[i+1]
                # Use networkx to find the shortest path length in meters
                length = nx.shortest_path_length(G, source=source_node, target=target_node, weight='length')
                total_length += length
            
            # Calculate walking time in minutes
            time_minutes = total_length / walking_speed_mpm
            results[route_id] = {'length': total_length, 'time': time_minutes}

            print(f"Route {route_id}:")
            print(f"  Path Length: {total_length:.0f} meters")
            print(f"  Calculation: {total_length:.0f} meters / {walking_speed_mpm} m/min")
            print(f"  Estimated Time: {time_minutes:.1f} minutes\n")

        except Exception:
            # Assign a very high value to any route that fails calculation
            results[route_id] = {'length': float('inf'), 'time': float('inf')}
    
    # Determine the fastest route by finding the minimum time
    if not any(v['time'] != float('inf') for v in results.values()):
        print("No routes could be calculated.", file=sys.stderr)
        return
        
    fastest_route_id = min(results, key=lambda k: results[k]['time'])
    fastest_time = results[fastest_route_id]['time']

    print("---")
    print("Conclusion:")
    print(f"The fastest route is Route {fastest_route_id} with an estimated walking time of {fastest_time:.1f} minutes.")


if __name__ == "__main__":
    solve_walking_route()