import time
from geopy.geocoders import Nominatim
from geopy.distance import geodesic
import sys

def get_coordinates(geolocator, address, cache):
    """
    Get coordinates for a given address, using a cache to avoid repeated requests.
    Includes a delay to respect API usage policies.
    """
    if address in cache:
        return cache[address]
    
    try:
        # Nominatim requires a small delay between requests
        time.sleep(1) 
        location = geolocator.geocode(address)
        if location:
            coords = (location.latitude, location.longitude)
            cache[address] = coords
            return coords
        else:
            print(f"Warning: Could not geocode '{address}'. Skipping.", file=sys.stderr)
            return None
    except Exception as e:
        print(f"Error geocoding '{address}': {e}", file=sys.stderr)
        return None

def calculate_total_distance(geolocator, waypoints, cache):
    """
    Calculate the total geodesic distance for a list of waypoints.
    """
    total_distance = 0
    # Convert waypoint names to coordinates
    coords_list = [get_coordinates(geolocator, wp, cache) for wp in waypoints]
    
    # Filter out any waypoints that couldn't be geocoded
    valid_coords = [c for c in coords_list if c is not None]
    
    if len(valid_coords) < 2:
        return float('inf') # Not a valid path

    # Calculate distance between consecutive points
    for i in range(len(valid_coords) - 1):
        point1 = valid_coords[i]
        point2 = valid_coords[i+1]
        total_distance += geodesic(point1, point2).meters
        
    return total_distance

def main():
    """
    Main function to find the fastest walking route.
    """
    print("Initializing geolocator and calculating route distances...")
    # Initialize geolocator with a unique user_agent
    geolocator = Nominatim(user_agent="route_comparison_tool_1.0")
    
    # A cache to store coordinates and avoid re-querying the same location
    coord_cache = {}

    start_point = "Guildhall, London"
    end_point = "St Paul's Cathedral, London"

    # Define the key waypoints for each route based on the description
    routes = {
        'A': [start_point, "Foster Lane and Gresham Street, London", end_point],
        'B': [start_point, "Rotunda, London Wall, London", "St. Martin's Le Grand, London", end_point],
        'C': [start_point, "Princes Street and Lothbury, London", "Queen Victoria Street, London", end_point],
        'D': [start_point, "St Swithin's Lane, London", "Cannon Street Station, London", "High Timber Street, London", end_point],
        'E': [start_point, "Wood Street and St Giles Terrace, London", "Aldersgate Street, London", end_point]
    }

    results = {}
    print("\nEstimated Route Distances:")
    for route_id, waypoints in routes.items():
        distance = calculate_total_distance(geolocator, waypoints, coord_cache)
        results[route_id] = distance
        print(f"Route {route_id}: {distance:.0f} meters")

    if not results:
        print("Could not calculate any route distances.")
        return

    # Find the route with the minimum distance
    fastest_route = min(results, key=results.get)
    min_distance = results[fastest_route]
    
    print("\n---")
    print(f"The fastest route is {fastest_route} with an estimated distance of {min_distance:.0f} meters.")
    print("---\n")


if __name__ == '__main__':
    main()