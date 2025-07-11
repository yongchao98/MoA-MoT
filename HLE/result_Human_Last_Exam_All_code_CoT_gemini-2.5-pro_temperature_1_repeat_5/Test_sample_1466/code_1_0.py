import math

def find_fastest_route():
    """
    Analyzes walking routes from Guildhall to St Paul's Cathedral,
    considering a road closure on Cheapside, to find the fastest option.
    """
    # Average walking speed is ~5 km/h, which is 5000 meters / 60 minutes.
    WALKING_SPEED_METERS_PER_MINUTE = 5000 / 60

    # Route data based on map analysis. Distances are estimated for each valid detour.
    # The closure on Cheapside is between Grocers' Hall Court and Gutter Lane.
    # Route A enters Cheapside from Foster Lane, which is inside the closed section.
    routes = {
        'A': {
            'description': "Impossible route. It enters Cheapside within the closed-off section.",
            'is_possible': False,
            'distance_m': None
        },
        'B': {
            'description': "Northern detour via London Wall and Aldersgate St.",
            'is_possible': True,
            'distance_m': 1150
        },
        'C': {
            'description': "Long southern detour via Queen Victoria St.",
            'is_possible': True,
            'distance_m': 1500
        },
        'D': {
            'description': "Extremely long southern detour via Upper Thames St.",
            'is_possible': True,
            'distance_m': 2200
        },
        'E': {
            'description': "Direct northern detour via Barbican and Aldersgate St.",
            'is_possible': True,
            'distance_m': 1050
        }
    }

    fastest_route_id = None
    min_time = float('inf')

    print("Analyzing walking routes from Guildhall to St Paul's Cathedral...")
    print(f"Closure: Cheapside between Grocers' Hall Court and Gutter Lane.\n")
    print("-" * 40)

    for route_id, data in routes.items():
        print(f"Route {route_id}:")
        print(f"  Description: {data['description']}")
        if data['is_possible']:
            distance = data['distance_m']
            # Time (min) = Distance (m) / Speed (m/min)
            time_minutes = distance / WALKING_SPEED_METERS_PER_MINUTE
            
            # Print the equation as requested
            print(f"  Calculation: {distance} meters / {WALKING_SPEED_METERS_PER_MINUTE:.1f} m/min = {time_minutes:.1f} minutes")

            if time_minutes < min_time:
                min_time = time_minutes
                fastest_route_id = route_id
        else:
            print("  Status: Invalid due to road closure.")
        print("-" * 40)

    print("\n[CONCLUSION]")
    if fastest_route_id:
        print(f"The fastest plausible route is '{fastest_route_id}'.")
        print(f"It has an estimated walking time of {min_time:.1f} minutes.")
    else:
        print("No possible routes were found.")

if __name__ == '__main__':
    find_fastest_route()