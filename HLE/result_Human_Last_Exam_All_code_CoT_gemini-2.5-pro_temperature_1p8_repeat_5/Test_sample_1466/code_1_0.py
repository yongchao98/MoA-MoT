import math

def find_fastest_route():
    """
    Analyzes walking routes from Guildhall to St Paul's Cathedral,
    calculates walking time for each, and determines the fastest option.
    """
    # Average walking speed is approximately 5 km/h.
    # 5 km/h = 5000 meters / 60 minutes
    walking_speed_meters_per_minute = 5000 / 60

    # Route data: distance estimations in meters based on map analysis.
    # A is invalid, D is suboptimal.
    routes = {
        'A': {'distance_m': -1, 'valid': False, 'notes': 'Invalid route, uses the closed section of Cheapside.'},
        'B': {'distance_m': 1100, 'valid': True, 'notes': 'Westward detour via Aldersgate St.'},
        'C': {'distance_m': 1200, 'valid': True, 'notes': 'Southern detour via Queen Victoria St.'},
        'D': {'distance_m': 1800, 'valid': False, 'notes': 'Suboptimal route, excessively long and complex.'},
        'E': {'distance_m': 1150, 'valid': True, 'notes': 'North-west detour via Barbican, slightly less direct than B.'}
    }

    fastest_route = None
    min_time = float('inf')

    print("Analyzing walking routes from Guildhall to St Paul's Cathedral...")
    print(f"Assumed walking speed: {round(walking_speed_meters_per_minute, 1)} meters/minute.")
    print("-" * 60)

    # Process routes in alphabetical order for consistent output
    for route_name in sorted(routes.keys()):
        route = routes[route_name]
        print(f"Route {route_name}: {route['notes']}")

        if route['valid']:
            distance = route['distance_m']
            # Calculate time and display the equation
            time_minutes = distance / walking_speed_meters_per_minute
            print(f"  Calculation: {distance} meters / {round(walking_speed_meters_per_minute, 1)} m/min = {round(time_minutes, 1)} minutes")

            if time_minutes < min_time:
                min_time = time_minutes
                fastest_route = route_name
        print("-" * 60)

    print("\n[Conclusion]")
    print(f"Comparing the valid routes, Route {fastest_route} is the fastest option.")

if __name__ == "__main__":
    find_fastest_route()