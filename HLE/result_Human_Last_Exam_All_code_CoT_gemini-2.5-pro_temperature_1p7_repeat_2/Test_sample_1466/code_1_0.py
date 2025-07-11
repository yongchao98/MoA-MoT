import requests
import json

def calculate_fastest_route():
    """
    This script calculates the fastest walking route from a set of options
    by querying a public routing API for an estimated travel time.
    """
    
    # Step 1: Define waypoints (longitude, latitude) for each route.
    # These points represent the key turns and locations in each route description.
    # The more detailed the waypoints, the more accurate the route calculation.
    waypoints = {
        'A': [
            (-0.0919, 51.5155), # Guildhall
            (-0.0924, 51.5149), # on Gresham St
            (-0.0954, 51.5147), # Gresham St/Foster Ln
            (-0.0957, 51.5144), # Foster Ln/Cheapside (This is west of the Gutter Ln closure)
            (-0.0984, 51.5138)  # St Paul's Cathedral
        ],
        'B': [
            (-0.0919, 51.5155), # Guildhall
            (-0.0905, 51.5173), # Coleman St/London Wall (Long detour north)
            (-0.0963, 51.5178), # Rotunda/Aldersgate St
            (-0.0968, 51.5149), # St Martin's Le Grand
            (-0.0984, 51.5138)  # St Paul's Cathedral
        ],
        'C': [
            (-0.0919, 51.5155), # Guildhall
            (-0.0894, 51.5142), # Lothbury/Princes St (Long detour south-east)
            (-0.0899, 51.5130), # Mansion House St
            (-0.0955, 51.5118), # Queen Victoria St/Cannon St
            (-0.0984, 51.5138)  # St Paul's Cathedral
        ],
        'D': [ 
            (-0.0919, 51.5155), # Guildhall
            (-0.0883, 51.5123), # Lombard St/Cannon St (Very long detour south)
            (-0.0909, 51.5101), # Upper Thames St
            (-0.0950, 51.5112), # near Queen Victoria St
            (-0.0984, 51.5138)  # St Paul's Cathedral
        ],
         'E': [ 
            (-0.0919, 51.5155), # Guildhall
            (-0.0934, 51.5181), # Wood St/Barbican Highwalk (Long detour north)
            (-0.0963, 51.5178), # Aldersgate St
            (-0.0984, 51.5138)  # St Paul's Cathedral
        ]
    }

    # Step 2: Define a function to call the OSRM routing API.
    def get_walking_duration_minutes(route_waypoints):
        """Fetches walking duration for a list of waypoints from OSRM API."""
        coordinates = ";".join([f"{lon},{lat}" for lon, lat in route_waypoints])
        url = f"http://router.project-osrm.org/route/v1/foot/{coordinates}?overview=false"
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            if data.get('code') == 'Ok':
                duration_seconds = data['routes'][0]['duration']
                return duration_seconds / 60
            return float('inf')
        except (requests.RequestException, KeyError, IndexError):
            return float('inf')

    # Step 3: Calculate duration for each route and store it.
    durations = {}
    print("Calculating walking times...")
    for route_id, route_waypoints in sorted(waypoints.items()):
        durations[route_id] = get_walking_duration_minutes(route_waypoints)
    
    # Step 4: Find the minimum duration and the corresponding route.
    # We filter out any routes that failed calculation (returned infinity).
    valid_durations = {k: v for k, v in durations.items() if v != float('inf')}
    if not valid_durations:
        print("Could not calculate walking times for any route.")
        return

    fastest_route_id = min(valid_durations, key=valid_durations.get)
    
    # Step 5: Print the results in a clear format.
    print("\n--- Route Comparison ---")
    for route_id, time_min in sorted(durations.items()):
        if time_min == float('inf'):
            print(f"Route {route_id}: Could not be calculated.")
        else:
            # Output the duration for each route as part of the comparison
            print(f"Route {route_id}: {time_min:.1f} minutes")
    
    print("\n--- Conclusion ---")
    print(f"The fastest route is '{fastest_route_id}' as it has the shortest walking time.")


if __name__ == '__main__':
    calculate_fastest_route()