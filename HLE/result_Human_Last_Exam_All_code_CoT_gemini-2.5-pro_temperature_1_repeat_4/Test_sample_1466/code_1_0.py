import openrouteservice
from geopy.geocoders import Nominatim
import time
import sys

def calculate_fastest_route():
    """
    Calculates the walking time for several predefined routes from Guildhall to St Paul's Cathedral
    and determines the fastest one.
    """
    # IMPORTANT: You must get a free API key from https://openrouteservice.org/
    # and replace "YOUR_API_KEY" with your actual key.
    api_key = 'YOUR_API_KEY' 
    
    if api_key == 'YOUR_API_KEY' or not api_key:
        print("Execution stopped: Please replace 'YOUR_API_KEY' in the script with your actual openrouteservice API key.", file=sys.stderr)
        return

    try:
        client = openrouteservice.Client(key=api_key)
        geolocator = Nominatim(user_agent="london_route_calculator_project")
    except Exception as e:
        print(f"Failed to initialize services. Check API key and network. Error: {e}", file=sys.stderr)
        return

    # Define key waypoints for each route to capture its general path.
    # These waypoints represent the essence of the described detours.
    route_waypoints = {
        'A': ["Guildhall, London", "Gresham St at Foster Ln, London", "St Paul's Cathedral, London"],
        'B': ["Guildhall, London", "London Wall at Coleman Street, London", "Aldersgate St at London Wall, London", "St Paul's Cathedral, London"],
        'C': ["Guildhall, London", "Lothbury at Princes St, London", "Queen Victoria Street at Cannon Street, London", "St Paul's Cathedral, London"],
        'D': ["Guildhall, London", "Lombard St, London", "Upper Thames Street, London", "St Paul's Cathedral, London"],
        'E': ["Guildhall, London", "Wood St near London Wall, London", "St. Martin's Le Grand, London", "St Paul's Cathedral, London"]
    }

    results = {}
    print("Calculating walking times for all routes...")
    print("-" * 40)

    for route_id, waypoints in route_waypoints.items():
        print(f"Processing Route {route_id}...")
        coords_list = []
        valid_waypoints = True
        for wp in waypoints:
            try:
                # Add a delay between geocoding requests to respect usage policies
                time.sleep(1) 
                location = geolocator.geocode(wp, timeout=10)
                if location:
                    coords_list.append((location.longitude, location.latitude))
                else:
                    print(f"Warning: Could not find coordinates for '{wp}'.", file=sys.stderr)
                    valid_waypoints = False
                    break
            except Exception as e:
                print(f"Error geocoding '{wp}': {e}", file=sys.stderr)
                valid_waypoints = False
                break
        
        if not valid_waypoints:
            print(f"Skipping Route {route_id} due to waypoint error.\n" + "-" * 40)
            continue

        try:
            # Get the route from the ORS API
            route_info = client.directions(
                coordinates=coords_list,
                profile='foot-walking'
            )
            # Duration is in seconds, convert to minutes
            duration_minutes = route_info['routes'][0]['summary']['duration'] / 60
            results[route_id] = duration_minutes
            # Output the number for this route
            print(f"Route {route_id} calculated time: {duration_minutes:.1f} minutes")

        except Exception as e:
            print(f"Could not calculate Route {route_id}. Error: {e}", file=sys.stderr)
        
        print("-" * 40)

    if not results:
        print("Could not calculate any routes. Please check your API key and waypoint names.", file=sys.stderr)
        return

    # Find the fastest route
    fastest_route_id = min(results, key=results.get)
    fastest_time = results[fastest_route_id]

    print("\n--- Final Results ---")
    print("The estimated walking time for each route is:")
    # Output each number in the final comparison
    for route_id, duration in sorted(results.items()):
        print(f"Route {route_id}: {duration:.1f} minutes")

    print(f"\nThe fastest walking route is Route {fastest_route_id} with an estimated time of {fastest_time:.1f} minutes.")
    print(f"<<<{fastest_route_id}>>>")


if __name__ == '__main__':
    calculate_fastest_route()