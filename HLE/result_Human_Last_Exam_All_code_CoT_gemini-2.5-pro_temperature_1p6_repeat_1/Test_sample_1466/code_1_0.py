import osmnx as ox
import networkx as nx

def solve():
    """
    Calculates the fastest walking route from Guildhall to St Paul's Cathedral
    by evaluating five given options.
    """
    # Use a cached graph for speed if available, otherwise fetch and cache
    try:
        G = ox.load_graphml("london_walk.graphml")
    except FileNotFoundError:
        print("Downloading London street map data. This may take a moment...")
        G = ox.graph_from_place('City of London, UK', network_type='walk')
        ox.save_graphml(G, "london_walk.graphml")

    # Average walking speed: 5 km/h = 5000 m/60 min ~= 83.3 m/min
    WALK_SPEED_MPM = 5000 / 60

    routes = {
        'A': [
            (51.5156, -0.0919), # Start: Guildhall
            (51.5152, -0.0935), # On Gresham St
            (51.5149, -0.0951), # On Foster Ln
            (51.5140, -0.0975), # On New Change
            (51.5138, -0.0984)  # End: St Paul's
        ],
        'B': [
            (51.5156, -0.0919), # Start: Guildhall
            (51.5173, -0.0908), # On London Wall
            (51.5181, -0.0963), # On Aldersgate St
            (51.5142, -0.0972), # On St Martin's Le Grand
            (51.5138, -0.0984)  # End: St Paul's
        ],
        'C': [
            (51.5156, -0.0919), # Start: Guildhall
            (51.5146, -0.0898), # On Lothbury
            (51.5130, -0.0908), # On Queen Victoria St
            (51.5123, -0.0945), # Near Cannon St
            (51.5138, -0.0984)  # End: St Paul's
        ],
        'D': [
            (51.5156, -0.0919), # Start: Guildhall
            (51.5146, -0.0898), # On Lothbury
            (51.5126, -0.0886), # On Cannon St
            (51.5115, -0.0914), # Near Dowgate Hill
            (51.5114, -0.0942), # On Upper Thames St
            (51.5131, -0.0968), # On Friday St
            (51.5138, -0.0984)  # End: St Paul's
        ],
        'E': [
            (51.5156, -0.0919), # Start: Guildhall
            (51.5164, -0.0946), # On Wood St
            (51.5179, -0.0950), # On Barbican Highwalk
            (51.5181, -0.0963), # On Aldersgate St
            (51.5142, -0.0972), # Near Cheapside/New Change
            (51.5138, -0.0984)  # End: St Paul's
        ]
    }

    results = {}

    print("Calculating walking times for all routes...\n")
    for name, waypoints_coords in routes.items():
        try:
            # Find the graph nodes closest to our waypoint coordinates
            waypoint_nodes = ox.nearest_nodes(G, [p[1] for p in waypoints_coords], [p[0] for p in waypoints_coords])

            total_length = 0
            segment_lengths = []
            # Calculate path length between each consecutive waypoint
            for i in range(len(waypoint_nodes) - 1):
                length = nx.shortest_path_length(G, waypoint_nodes[i], waypoint_nodes[i+1], weight='length')
                segment_lengths.append(length)
                total_length += length
            
            # Calculate total time
            total_time = total_length / WALK_SPEED_MPM
            results[name] = total_time
            
            # Format the equation string
            segments_str = " + ".join([f"{l:.0f}m" for l in segment_lengths])
            print(f"Route {name}:")
            print(f"Equation: {segments_str} = {total_length:.0f} meters")
            print(f"Estimated Time: {total_time:.1f} minutes\n")

        except nx.NetworkXNoPath:
            print(f"Could not calculate Route {name} as a valid path was not found.\n")
            results[name] = float('inf')


    if not results:
        print("Could not calculate any routes.")
        return

    # Find the fastest route
    fastest_route = min(results, key=results.get)
    fastest_time = results[fastest_route]
    
    print("-" * 30)
    print(f"The fastest route is Route {fastest_route} with an estimated time of {fastest_time:.1f} minutes.")
    print(f"<<<{fastest_route}>>>")

solve()