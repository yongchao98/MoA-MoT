import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance in meters between two points
    on the earth (specified in decimal degrees).
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.asin(math.sqrt(a))
    
    # Radius of earth in meters.
    r_meters = 6371000
    return c * r_meters

def solve_walking_route():
    """
    Analyzes the walking routes, calculates their distances and times,
    and determines the fastest option.
    """
    print("To find the fastest walking route while Cheapside is closed, I will model the given routes using geographic coordinates and calculate their lengths.")
    print("The fastest route will be the one with the shortest total distance. I'll use the Haversine formula to calculate distances between points and assume a constant walking speed.\n")
    
    # Step 1: Define the coordinates for key points in each route.
    waypoints = {
        'Guildhall': (51.5155, -0.0920),
        'St_Pauls': (51.5138, -0.0984),
        # Route B specific waypoints
        'Coleman_LondonWall': (51.5175, -0.0917),
        'Rotunda': (51.5173, -0.0957), # Junction of London Wall & Aldersgate St
        'StMartins_Cheapside': (51.5144, -0.0975),
        # Route C specific waypoints
        'Bank_Junction': (51.5134, -0.0886),
        'QV_Cannon': (51.5126, -0.0941), # Junction of Queen Victoria St & Cannon St
        # Route E specific waypoints
        'Wood_St_Junction': (51.5165, -0.0945),
    }

    # Define the routes as a sequence of waypoints and descriptions
    routes = {
        "A": {
            "path": [],
            "description": "Start at Guildhall -> Gresham St -> Foster Ln -> Cheapside -> ...",
            "comment": "Invalid. This route attempts to use the closed section of Cheapside."
        },
        "B": {
            "path": ['Guildhall', 'Coleman_LondonWall', 'Rotunda', 'StMartins_Cheapside', 'St_Pauls'],
            "description": "Start at Guildhall -> ... -> London Wall -> Aldersgate St -> ...",
            "comment": "A plausible northern detour."
        },
        "C": {
            "path": ['Guildhall', 'Bank_Junction', 'QV_Cannon', 'St_Pauls'],
            "description": "Start at Guildhall -> ... -> Princes St -> Queen Victoria St -> ...",
            "comment": "A plausible southern detour."
        },
        "D": {
            "path": [],
            "description": "Start at Guildhall -> ... -> Upper Thames St -> ...",
            "comment": "Invalid. This route is clearly circuitous and significantly longer."
        },
        "E": {
            "path": ['Guildhall', 'Wood_St_Junction', 'Rotunda', 'StMartins_Cheapside', 'St_Pauls'],
            "description": "Start at Guildhall -> ... -> Wood St -> Barbican -> Aldersgate St -> ...",
            "comment": "A more direct northern detour."
        }
    }
    
    # Step 2: Calculate distance and time for each valid route.
    # Average walking speed: 5 km/h = 83.33 meters/minute
    walking_speed_mpm = 83.33
    results = {}

    print("Step 2: Calculate the distance and estimated time for each valid route.\n")

    for route_name, route_data in routes.items():
        print(f"--- Analyzing Route {route_name} ---")
        print(f"Description: {route_data['description']}")
        print(f"Verdict: {route_data['comment']}")

        if not route_data["path"]:
            print("\n")
            continue

        total_distance = 0
        path = route_data["path"]
        segment_distances = []

        print("Calculation:")
        for i in range(len(path) - 1):
            p1_name, p2_name = path[i], path[i+1]
            p1_coords, p2_coords = waypoints[p1_name], waypoints[p2_name]
            
            segment_distance = haversine(p1_coords[0], p1_coords[1], p2_coords[0], p2_coords[1])
            segment_distances.append(segment_distance)
            total_distance += segment_distance
            
            print(f"  Distance({p1_name} -> {p2_name}) = {segment_distance:.0f} meters")

        walking_time = total_distance / walking_speed_mpm
        results[route_name] = {'distance': total_distance, 'time': walking_time}
        
        distance_eq = " + ".join([f"{d:.0f}" for d in segment_distances])
        print(f"Total Distance = {distance_eq} = {total_distance:.0f} meters")
        print(f"Estimated Time ({total_distance:.0f}m / {walking_speed_mpm:.2f}m/min) = {walking_time:.1f} minutes\n")

    # Step 3: Compare results and find the fastest route.
    if not results:
         print("No valid routes could be calculated.")
         return

    fastest_route_name = min(results, key=lambda k: results[k]['time'])
    fastest_route_info = results[fastest_route_name]

    print("--- Conclusion ---")
    print(f"Comparing the estimated times:")
    for name, info in sorted(results.items(), key=lambda item: item[1]['time']):
        print(f"  Route {name}: {info['time']:.1f} minutes")

    print(f"\nThe fastest route is Route {fastest_route_name} with an estimated time of {fastest_route_info['time']:.1f} minutes.")

solve_walking_route()
<<<E>>>