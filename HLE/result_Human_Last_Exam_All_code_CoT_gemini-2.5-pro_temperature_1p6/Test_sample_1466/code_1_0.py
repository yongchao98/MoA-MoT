import osmnx as ox
import networkx as nx

def find_fastest_route():
    """
    Calculates the walking time for five different routes from Guildhall to St Paul's Cathedral
    and determines which is the fastest.
    """
    print("Starting the calculation for the fastest walking route...\n")

    # Step 1: Define start and end coordinates
    guildhall_coord = (51.5155, -0.0917)
    st_pauls_coord = (51.5138, -0.0984)

    # Step 2: Define waypoints to represent each route's path
    routes_waypoints = {
        'A': [guildhall_coord, (51.5146, -0.0945), st_pauls_coord],  # Via Gresham St/Foster Ln
        'B': [guildhall_coord, (51.5173, -0.0954), (51.5152, -0.0991), st_pauls_coord], # Via London Wall/Aldersgate St
        'C': [guildhall_coord, (51.5135, -0.0890), st_pauls_coord], # Via Bank Junction
        'D': [guildhall_coord, (51.5125, -0.0875), (51.5098, -0.0947), st_pauls_coord], # Via Cannon St/Upper Thames
        'E': [guildhall_coord, (51.5181, -0.0927), (51.5173, -0.0954), st_pauls_coord]  # Via Barbican/Aldersgate St
    }

    try:
        # Step 3: Download the walking network graph for the area
        print("Downloading London street map data. This may take a moment...")
        # Define a bounding box large enough to contain all routes
        north, south, east, west = 51.521, 51.508, -0.085, -0.101
        G = ox.graph_from_bbox(north, south, east, west, network_type='walk')

        # Use an average walking speed of 4.5 km/h to calculate travel times
        G = ox.add_edge_speeds(G, hwy_speeds=4.5)
        G = ox.add_edge_travel_times(G)
        print("Map data processed.\n")

        # Step 4: Calculate the travel time for each route
        route_times = {}
        print("Calculating travel time for each route:")
        for name, waypoints in routes_waypoints.items():
            total_time_seconds = 0
            try:
                # Get the graph nodes closest to our waypoint coordinates
                nodes = ox.nearest_nodes(G, [c[1] for c in waypoints], [c[0] for c in waypoints])
                # Calculate the time for the path segments between consecutive waypoints
                for i in range(len(nodes) - 1):
                    path_time = nx.shortest_path_length(G, source=nodes[i], target=nodes[i+1], weight='travel_time')
                    total_time_seconds += path_time
                route_times[name] = total_time_seconds / 60  # Convert to minutes
            except nx.NetworkXNoPath:
                route_times[name] = float('inf') # Should not happen in a connected city graph

        # Step 5: Print the results and find the minimum
        for name, time in sorted(route_times.items()):
             print(f"Route {name} estimated time: {time:.2f} minutes")

        fastest_route_name = min(route_times, key=route_times.get)
        fastest_time = route_times[fastest_route_name]

        print(f"\nBased on the calculations, the fastest route is Route {fastest_route_name} with an estimated time of {fastest_time:.2f} minutes.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("This may be due to a network issue or missing libraries (e.g., osmnx, geopandas).")
        print("Please ensure you have an internet connection and have installed the required libraries with 'pip install osmnx'.")

if __name__ == '__main__':
    find_fastest_route()