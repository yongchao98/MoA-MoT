import osmnx as ox
import networkx as nx

def find_fastest_route():
    """
    This function calculates the walking time for several predefined routes
    from Guildhall to St Paul's Cathedral in London and identifies the fastest one.
    """
    print("Initializing walking route analysis...")
    print("This may take a moment as map data is downloaded.")

    try:
        # Average walking speed: 5 km/h, converted to meters per minute.
        # 5000 meters / 60 minutes
        walking_speed_mpm = 5000 / 60

        # Download the walking network for the area of interest in London.
        # Using a bounding box to ensure all routes are covered.
        north, south, east, west = 51.519, 51.510, -0.088, -0.100
        G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
        print("Street network for London downloaded successfully.\n")

        # Define waypoints (latitude, longitude) for each route based on the descriptions.
        start_point = (51.5155, -0.0919) # Guildhall
        end_point = (51.5138, -0.0984)   # St Paul's Cathedral

        routes = {
            "A": [start_point, (51.5151, -0.0960), (51.5141, -0.0980), end_point],
            "B": [start_point, (51.5176, -0.0955), (51.5144, -0.0976), end_point],
            "C": [start_point, (51.5131, -0.0905), (51.5126, -0.0945), end_point],
            "D": [start_point, (51.5132, -0.0919), (51.5114, -0.0935), (51.5125, -0.0973), end_point],
            "E": [start_point, (51.5173, -0.0945), (51.5176, -0.0955), (51.5144, -0.0976), end_point],
        }

        route_times = {}

        print("--- Calculating Route Durations ---")
        for name, waypoints in routes.items():
            total_distance = 0
            # Get the graph nodes closest to the route's waypoints.
            try:
                nodes = ox.distance.nearest_nodes(G, [p[1] for p in waypoints], [p[0] for p in waypoints])

                # Calculate the shortest path distance between each pair of consecutive waypoints.
                for i in range(len(nodes) - 1):
                    path_length = nx.shortest_path_length(G, source=nodes[i], target=nodes[i+1], weight='length')
                    total_distance += path_length
                
                # Calculate time and print the equation.
                time_in_minutes = total_distance / walking_speed_mpm
                route_times[name] = time_in_minutes
                
                print(f"Route {name}: {total_distance:.0f} meters / {walking_speed_mpm:.2f} m/min = {time_in_minutes:.1f} minutes")

            except nx.NetworkXNoPath:
                print(f"Route {name}: Could not find a valid path.")
                route_times[name] = float('inf')
        
        # Find the route with the minimum time.
        if not route_times:
            print("Could not calculate times for any route.")
            return

        fastest_route_name = min(route_times, key=route_times.get)
        fastest_time = route_times[fastest_route_name]

        print("\n--- Conclusion ---")
        print(f"The fastest route is option {fastest_route_name}, with an estimated walking time of {fastest_time:.1f} minutes.")
        
        # Print the final answer in the required format.
        print(f"\n<<<{fastest_route_name}>>>")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and the required libraries (osmnx, networkx) are installed (`pip install osmnx`).")

if __name__ == "__main__":
    find_fastest_route()