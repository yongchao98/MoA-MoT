# To run this code, you may need to install the 'geopy' library.
# You can do this by running the following command in your terminal:
# pip install geopy

import geopy.distance

def solve_and_print_route():
    """
    This script calculates the fastest walking route from Guildhall to St Paul's Cathedral
    with a road closure on Cheapside. It analyzes several given routes by calculating
    their total length based on GPS coordinates.
    """
    
    # Step 1: Define GPS coordinates for key points along each suggested route.
    waypoints = {
        'Start: Guildhall': (51.5155, -0.0919),
        'End: St Pauls Cathedral': (51.5138, -0.0984),

        # Route B Points
        'Walk south on Basinghall St': (51.5151, -0.0925),
        'Turn left onto Masons Ave': (51.5153, -0.0905),
        'Turn left onto Coleman St': (51.5156, -0.0900),
        'Turn left onto London Wall': (51.5165, -0.0910),
        'To Aldersgate Rotunda': (51.5169, -0.0950),
        'Continue onto St Martins Le Grand': (51.5152, -0.0967),
        'To New Change': (51.5142, -0.0975),

        # Route C Points
        'Gresham St to Lothbury': (51.5147, -0.0895),
        'Lothbury to Princes St': (51.5137, -0.0886),
        'To Mansion House St': (51.5131, -0.0894),
        'Onto Queen Victoria St': (51.5124, -0.0917),
        'Turn left to stay on Queen Victoria St': (51.5120, -0.0933),
        'To Cannon St': (51.5117, -0.0950),
        
        # Route D Points
        'Turn left onto Lombard St': (51.5124, -0.0864),
        'Turn right onto St Swithins Ln': (51.5119, -0.0885),
        'Turn right onto Cannon St': (51.5116, -0.0893),
        'Turn left onto Dowgate Hill': (51.5109, -0.0912),
        'Turn right towards Upper Thames St': (51.5106, -0.0915),
        'Turn right onto Upper Thames St': (51.5102, -0.0928),
        'Turn left onto High Timber St': (51.5102, -0.0953),
        'Turn right towards Fye Foot Ln':(51.5118, -0.0963),
        'Turn left onto Friday St': (51.5126, -0.0970),

        # Route E Points
        'Walk north on Basinghall St': (51.5158, -0.0919),
        'Turn right towards Wood St': (51.5173, -0.0934),
        'Turn right onto Barbican Highwalk': (51.5178, -0.0945),
        'Turn right onto Aldersgate St': (51.5173, -0.0958),
    }

    # Step 2: Define the sequence of waypoints for each route.
    routes = {
        'A': 'Invalid (uses closed section of Cheapside)',
        'B': ['Start: Guildhall', 'Walk south on Basinghall St', 'Turn left onto Masons Ave', 'Turn left onto Coleman St', 'Turn left onto London Wall', 'To Aldersgate Rotunda', 'Continue onto St Martins Le Grand', 'To New Change', 'End: St Pauls Cathedral'],
        'C': ['Start: Guildhall', 'Gresham St to Lothbury', 'Lothbury to Princes St', 'To Mansion House St', 'Onto Queen Victoria St', 'Turn left to stay on Queen Victoria St', 'To Cannon St', 'End: St Pauls Cathedral'],
        'D': ['Start: Guildhall', 'Gresham St to Lothbury', 'Lothbury to Princes St', 'Turn left onto Lombard St', 'Turn right onto St Swithins Ln', 'Turn right onto Cannon St', 'Turn left onto Dowgate Hill', 'Turn right towards Upper Thames St', 'Turn right onto Upper Thames St', 'Turn left onto High Timber St', 'Turn right towards Fye Foot Ln', 'Turn left onto Friday St', 'End: St Pauls Cathedral'],
        'E': ['Start: Guildhall', 'Walk north on Basinghall St', 'Turn right towards Wood St', 'Turn right onto Barbican Highwalk', 'Turn right onto Aldersgate St', 'To Aldersgate Rotunda', 'Continue onto St Martins Le Grand', 'To New Change', 'End: St Pauls Cathedral']
    }

    # Step 3 & 4: Calculate distances for valid routes.
    print("--- Route Analysis ---")
    print("Route A is invalid as it passes through the road closure on Cheapside.")
    print("Calculating distances for the other routes...\n")
    
    results = {}
    for route_name, route_points in routes.items():
        if route_name == 'A':
            continue
        
        total_distance = 0
        for i in range(len(route_points) - 1):
            point1 = waypoints[route_points[i]]
            point2 = waypoints[route_points[i+1]]
            total_distance += geopy.distance.great_circle(point1, point2).km
        results[route_name] = total_distance

    print("--- Route Distances ---")
    sorted_results = sorted(results.items(), key=lambda item: item[1])
    for name, dist in sorted_results:
        print(f"Route {name}: {dist:.2f} km")
    
    best_route_name = sorted_results[0][0]
    
    print(f"\n--- Fastest Route Found: Route {best_route_name} ---")
    print("Below is the segment-by-segment distance calculation for the fastest route:")

    # Step 5 & 6: Print the detailed breakdown for the winning route as an "equation".
    best_route_points = routes[best_route_name]
    equation_parts = []
    
    for i in range(len(best_route_points) - 1):
        p1_name = best_route_points[i]
        p2_name = best_route_points[i+1]
        dist_segment = geopy.distance.great_circle(waypoints[p1_name], waypoints[p2_name]).km
        equation_parts.append(f"{dist_segment:.3f}")
        
    print(" + ".join(equation_parts) + f" = {results[best_route_name]:.3f} km")

    print(f"\nBased on the analysis, Route {best_route_name} is the shortest and therefore fastest path.")
    print(f"<<<{best_route_name}>>>")


solve_and_print_route()