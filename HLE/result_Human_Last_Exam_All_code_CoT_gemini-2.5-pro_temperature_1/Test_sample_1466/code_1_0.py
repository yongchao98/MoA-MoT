import osmnx as ox
import networkx as nx

def solve_walking_route():
    """
    Calculates the fastest walking route from Guildhall to St Paul's Cathedral
    given five route options and a road closure.

    Prerequisites: You need to install osmnx and its dependencies.
    'pip install osmnx'
    """
    try:
        # 1. Model the City: Download the walking network for the City of London.
        place_name = "City of London, London, UK"
        G = ox.graph_from_place(place_name, network_type='walk')
    except Exception as e:
        print(f"Could not download the map data. Please check your internet connection. Error: {e}")
        return

    # Define Start and End points
    try:
        start_loc = "Guildhall, London, UK"
        end_loc = "St Paul's Cathedral, London, UK"
        start_coords = ox.geocode(start_loc)
        end_coords = ox.geocode(end_loc)
    except Exception as e:
        print(f"Could not find start or end locations. Error: {e}")
        return

    # 2. Define the Routes: Specify key waypoints for each route option.
    # These points guide the path calculation to match the route descriptions.
    routes_waypoints = {
        'A': ["Gresham St and Foster Ln", "Cheapside and New Change"],
        'B': ["London Wall and Coleman St", "Aldersgate St and London Wall", "St Martin's Le Grand"],
        'C': ["Gresham St and Princes St", "Queen Victoria St and Cannon St"],
        'D': ["Lothbury and Princes St", "Cannon St and St Swithin's Ln", "Upper Thames St"],
        'E': ["Wood St and St Giles Terrace", "Aldersgate St and Barbican Highwalk", "St Martin's Le Grand"]
    }

    results = {}
    
    print("Calculating walking times for each route...\n")

    # 3. & 4. Calculate Distance and Time for each route
    # Average walking speed: 5 km/h = 5000 meters / 60 minutes = 83.33 meters per minute
    walking_speed_mpm = 5000 / 60

    for route_id, waypoints_text in routes_waypoints.items():
        try:
            # Geocode text waypoints to (latitude, longitude)
            waypoints_coords = [ox.geocode(point + ", London, UK") for point in waypoints_text]
            
            # Combine all points for the route calculation
            all_points = [start_coords] + waypoints_coords + [end_coords]
            
            total_length = 0
            # Calculate path length between consecutive points
            for i in range(len(all_points) - 1):
                # Find the nearest map nodes to our geocoded points
                origin_node = ox.distance.nearest_nodes(G, all_points[i][1], all_points[i][0])
                destination_node = ox.distance.nearest_nodes(G, all_points[i+1][1], all_points[i+1][0])
                
                # Calculate the shortest path length on the graph
                path_length = nx.shortest_path_length(G, origin_node, destination_node, weight='length')
                total_length += path_length
            
            # Convert distance to time
            time_minutes = total_length / walking_speed_mpm
            results[route_id] = {'distance': total_length, 'time': time_minutes}

        except Exception:
            # If a waypoint can't be found or a path doesn't exist, mark it as invalid.
            results[route_id] = {'distance': float('inf'), 'time': float('inf')}

    # 5. Compare and Conclude
    print("--- Route Analysis ---")
    for route_id, data in sorted(results.items()):
        if data['time'] != float('inf'):
            print(f"Route {route_id}: Distance = {data['distance']:.0f} m, Estimated Time = {data['time']:.1f} min")
        else:
            print(f"Route {route_id}: Could not calculate path.")
    
    # Find the fastest route
    if not results or all(data['time'] == float('inf') for data in results.values()):
        print("\nCould not determine the fastest route.")
        return
        
    fastest_route_id = min(results, key=lambda k: results[k]['time'])
    fastest_data = results[fastest_route_id]

    print("\n--- Conclusion ---")
    print(f"The fastest option is Route {fastest_route_id}.")
    
    # Display the final equation as requested
    print("\nFinal Time Calculation for Fastest Route:")
    distance = fastest_data['distance']
    time = fastest_data['time']
    print(f"Walking Speed = 5000 meters / 60 minutes = {walking_speed_mpm:.1f} m/min")
    print(f"Route {fastest_route_id} Distance = {distance:.0f} meters")
    print(f"Estimated Time = {distance:.0f} / {walking_speed_mpm:.1f} = {time:.1f} minutes")

if __name__ == '__main__':
    solve_walking_route()