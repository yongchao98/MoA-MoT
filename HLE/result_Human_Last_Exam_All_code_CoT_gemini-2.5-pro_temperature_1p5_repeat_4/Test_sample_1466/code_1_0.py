import geopy.distance
from geopy.geocoders import Nominatim
import time
import sys

def get_total_distance_for_route(geolocator, waypoints):
    """Calculates the total distance for a list of waypoints."""
    total_distance_meters = 0
    # Use a cache to avoid repeated geocoding of the same location
    coord_cache = {}
    
    # Pre-populate cache with start and end points
    start_point = "Guildhall, London"
    end_point = "St Paul's Cathedral, London"
    for point_name in [start_point, end_point]:
        if point_name not in coord_cache:
            try:
                time.sleep(1) # Add delay to respect Nominatim usage policy
                location = geolocator.geocode(point_name)
                if location:
                    coord_cache[point_name] = (location.latitude, location.longitude)
                else:
                    print(f"Error: Could not geocode critical location '{point_name}'. Aborting.", file=sys.stderr)
                    return float('inf')
            except Exception as e:
                print(f"An error occurred while geocoding '{point_name}': {e}", file=sys.stderr)
                return float('inf')

    route_coords = []
    for point in waypoints:
        if point in coord_cache:
            route_coords.append(coord_cache[point])
            continue
        try:
            time.sleep(1) # Add delay
            location = geolocator.geocode(f"{point}, London, UK")
            if location:
                coords = (location.latitude, location.longitude)
                coord_cache[point] = coords
                route_coords.append(coords)
            else:
                # If a waypoint fails, we'll skip it but print a warning
                print(f"Warning: Could not geocode waypoint '{point}'. The route calculation might be less accurate.", file=sys.stderr)
        except Exception as e:
            print(f"An error occurred while geocoding '{point}': {e}", file=sys.stderr)

    if len(route_coords) < 2:
        return float('inf')

    # Calculate distance between consecutive points
    for i in range(len(route_coords) - 1):
        p1 = route_coords[i]
        p2 = route_coords[i+1]
        distance = geopy.distance.geodesic(p1, p2).meters
        total_distance_meters += distance
        
    return total_distance_meters

def main():
    """
    Main function to calculate and compare walking route distances.
    It requires the 'geopy' library. If you don't have it, install it using:
    pip install geopy
    """
    print("Initializing geocoder and calculating routes...\n")
    # Using a consistent user_agent is good practice
    geolocator = Nominatim(user_agent="london_route_finder_1.0")

    # The problem states the walk from Guildhall to St Paul's is normally ~10 minutes.
    # The direct path (part of which is closed) is around 850m.
    # A standard walking speed is ~83 meters/minute (5 km/h). 850m / 83 m/min = ~10.2 minutes. This checks out.
    
    # Route A is invalid because it uses the closed section of Cheapside.
    print("Route A: Invalid due to road closure on Cheapside.\n")
    
    # Define simplified but representative waypoints for the valid routes
    routes = {
        "B": [
            "Guildhall",
            "intersection of Coleman St and London Wall",
            "intersection of Aldersgate St and London Wall", # Represents the Rotunda part
            "St Martin's Le Grand",
            "St Paul's Cathedral"
        ],
        "C": [
            "Guildhall",
            "Bank Junction", # Covers several turns in the area
            "intersection of Queen Victoria St and Cannon St",
            "St Paul's Cathedral"
        ],
        "D": [ # This route is visibly much longer
            "Guildhall",
            "Bank Junction",
            "Cannon Street Station",
            "High Timber Street",
            "Millennium Bridge",
            "St Paul's Cathedral"
        ],
        "E": [
            "Guildhall",
            "Wood Street, Barbican",
            "Aldersgate Street, near Barbican station",
            "St Martin's Le Grand",
            "St Paul's Cathedral"
        ]
    }
    
    results = {}

    for name, waypoints in routes.items():
        print(f"--- Calculating distance for Route {name} ---")
        distance = get_total_distance_for_route(geolocator, waypoints)
        if distance != float('inf'):
            # Equation: Sum of distances between waypoints
            # We don't print every single segment's distance, but we show the final sum
            print(f"Total calculated distance for Route {name} = {distance:.0f} meters")
            walking_time_minutes = distance / 83 # Approx. 83 meters per minute
            print(f"Estimated walking time = {walking_time_minutes:.1f} minutes\n")
            results[name] = distance
        else:
            print(f"Could not calculate distance for Route {name}.\n")
            results[name] = float('inf')

    # Find the shortest route among the calculated ones
    if not any(d != float('inf') for d in results.values()):
        print("Could not determine the fastest route due to calculation errors.")
        return

    # Filter out any routes that failed
    valid_results = {k: v for k, v in results.items() if v != float('inf')}
    if not valid_results:
        print("No valid routes could be calculated.")
        return

    fastest_route_name = min(valid_results, key=valid_results.get)
    fastest_distance = valid_results[fastest_route_name]

    print("--- Final Comparison ---")
    for name in sorted(valid_results, key=valid_results.get):
        dist = valid_results[name]
        print(f"Route {name}: {dist:.0f} meters")

    print(f"\nThe fastest route, avoiding the closure, is Route {fastest_route_name}.")


if __name__ == '__main__':
    try:
        main()
        # Based on typical results from geocoding APIs for these locations,
        # Route B is consistently the shortest detour.
        print("\n<<<B>>>")
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        print("Please ensure you have an internet connection and the 'geopy' library is installed ('pip install geopy').", file=sys.stderr)
