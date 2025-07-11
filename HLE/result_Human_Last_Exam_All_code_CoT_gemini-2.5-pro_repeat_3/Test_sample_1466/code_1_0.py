import osmnx as ox
import networkx as nx
import numpy as np

def solve_route_problem():
    """
    This function calculates the fastest walking route from a set of options
    by modeling the routes on a real-world street network from OpenStreetMap.
    """
    
    # Step 1: Define waypoint coordinates for each route.
    # Routes are defined by lists of (latitude, longitude) tuples from start to end.
    routes = {
        'A': [
            (51.5155, -0.0918),  # Guildhall
            (51.5152, -0.0957),  # Gresham St / Foster Ln
            (51.5148, -0.0960),  # Cheapside / Foster Ln
            (51.5141, -0.0975),  # Cheapside / New Change
            (51.5138, -0.0984)   # St Paul's Cathedral
        ],
        'B': [
            (51.5155, -0.0918),  # Guildhall
            (51.5173, -0.0917),  # London Wall / Coleman St
            (51.5179, -0.0960),  # Aldersgate St Roundabout
            (51.5146, -0.0970),  # St Martin's Le Grand / Cheapside
            (51.5138, -0.0984)   # St Paul's Cathedral
        ],
        'C': [
            (51.5155, -0.0918),  # Guildhall
            (51.5143, -0.0890),  # Lothbury / Princes St
            (51.5126, -0.0934),  # Queen Victoria St / Cannon St
            (51.5138, -0.0984)   # St Paul's Cathedral
        ],
        'D': [
            (51.5155, -0.0918),  # Guildhall
            (51.5123, -0.0883),  # Cannon St / St Swithin's Ln
            (51.5106, -0.0908),  # Upper Thames St / Dowgate Hill
            (51.5129, -0.0959),  # Queen Victoria St / Friday St
            (51.5138, -0.0984)   # St Paul's Cathedral
        ],
        'E': [
            (51.5155, -0.0918),  # Guildhall
            (51.5170, -0.0958),  # Wood St / St Giles Terrace
            (51.5186, -0.0970),  # Aldersgate St / Barbican Highwalk
            (51.5146, -0.0970),  # St Martin's Le Grand / Cheapside
            (51.5138, -0.0984)   # St Paul's Cathedral
        ]
    }

    # Step 2: Download the walking network graph for the area.
    print("Downloading street network for the City of London...")
    try:
        # Use a bounding box that comfortably contains all routes
        north, south, east, west = 51.520, 51.510, -0.085, -0.101
        G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
        
        # Add edge speeds and travel times for routing
        G = ox.add_edge_speeds(G)
        G = ox.add_edge_travel_times(G)
        print("Network download complete.\n")
    except Exception as e:
        print(f"Could not download map data. Please check your internet connection. Error: {e}")
        return

    # Step 3: Calculate the total walking time for each route.
    route_times = {}
    
    print("Calculating walking time for each route option:")
    for name, waypoints in routes.items():
        total_time_seconds = 0
        try:
            # Get the graph nodes closest to the waypoint coordinates
            lats = [p[0] for p in waypoints]
            lngs = [p[1] for p in waypoints]
            nodes = ox.nearest_nodes(G, X=lngs, Y=lats)
            
            # Calculate time between consecutive waypoints
            for i in range(len(nodes) - 1):
                # Use travel_time as the weight for the shortest path calculation
                path_time = nx.shortest_path_length(G, source=nodes[i], target=nodes[i+1], weight='travel_time')
                total_time_seconds += path_time
            
            route_times[name] = total_time_seconds
        except nx.NetworkXNoPath:
            # In case a path is not found between waypoints, assign a very high time
            route_times[name] = np.inf
            print(f"Could not calculate a complete path for Route {name}.")

    # Step 4: Display the results and identify the fastest route.
    print("\n--- Route Time Calculations ---")
    if not route_times:
        print("Could not calculate times for any route.")
        return
        
    for name, time_sec in sorted(route_times.items()):
        time_min = time_sec / 60
        # This output represents the "final equation" for each route
        print(f"Route {name} = {time_sec:.0f} seconds = {time_min:.2f} minutes")

    fastest_route_name = min(route_times, key=route_times.get)
    fastest_time_min = route_times[fastest_route_name] / 60
    
    print("\n--- Conclusion ---")
    print(f"The fastest route is {fastest_route_name} with an estimated walking time of {fastest_time_min:.2f} minutes.")
    
    # Final answer in the required format
    print(f"\n<<<{fastest_route_name}>>>")


if __name__ == '__main__':
    solve_route_problem()