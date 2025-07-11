# First, you may need to install the geopy library. You can do this by running:
# pip install geopy

import geopy.geocoders
from geopy.distance import geodesic
import time
import sys

def calculate_route_distance(geolocator, waypoints):
    """
    Geocodes a list of waypoints and calculates the total geodesic distance.
    Also prints the equation for the distance sum.
    """
    total_distance = 0.0
    segment_distances = []
    coords = []

    # Get coordinates for all waypoints
    for point in waypoints:
        try:
            # Adding a slight delay to respect the API's usage policy
            time.sleep(1)
            location = geolocator.geocode(point)
            if location:
                coords.append((location.latitude, location.longitude))
            else:
                print(f"    - Warning: Could not find location for '{point}'. Skipping route.", file=sys.stderr)
                return float('inf')
        except Exception as e:
            print(f"    - Error geocoding '{point}': {e}. Skipping route.", file=sys.stderr)
            return float('inf')

    # Calculate distance between each pair of consecutive waypoints
    for i in range(len(coords) - 1):
        segment_distance = geodesic(coords[i], coords[i+1]).meters
        segment_distances.append(segment_distance)
        total_distance += segment_distance
    
    # Print the equation representing the sum of segment distances
    distance_equation = " + ".join([f"{d:.0f}" for d in segment_distances])
    print(f"  Calculation: {distance_equation} = {total_distance:.0f} meters")

    return total_distance

# --- Main script ---
try:
    geolocator = geopy.geocoders.Nominatim(user_agent="london_route_finder_1.0")
except Exception as e:
    print(f"Failed to initialize geolocator: {e}", file=sys.stderr)
    sys.exit(1)


# Define the start and end points common to all routes
start_point = "Guildhall, London, UK"
end_point = "St Paul's Cathedral, London, UK"

# Define the waypoints that characterize each route
routes = {
    "A": [
        start_point,
        "intersection of Gresham St and Foster Ln, London",
        end_point
    ],
    "B": [
        start_point,
        "intersection of Coleman St and London Wall, A1211, London",
        "St Martin's Le Grand, near Aldersgate St, London",
        end_point
    ],
    "C": [
        start_point,
        "intersection of Gresham St and Princes St, London",
        "Queen Victoria Street, near Mansion House Station, London",
        end_point
    ],
    "D": [
        start_point,
        "intersection of Princes St and Threadneedle St, London",
        "Cannon Street, near St Swithin's Lane, London",
        end_point
    ],
    "E": [
        start_point,
        "Barbican Highwalk, near Wood Street, London",
        "Aldersgate Street, A1, near Barbican station, London",
        end_point
    ]
}

route_distances = {}
print("Calculating estimated walking distances for each route...\n")

for name, waypoints in routes.items():
    print(f"Analyzing Route {name}:")
    distance = calculate_route_distance(geolocator, waypoints)
    route_distances[name] = distance
    print(f"  Total estimated distance for Route {name}: {distance:.0f} meters\n")


# Find and announce the fastest route (shortest distance)
valid_routes = {k: v for k, v in route_distances.items() if v != float('inf')}

if not valid_routes:
    print("Could not calculate distances for any routes. Please check your connection or waypoint names.")
else:
    # Sort routes by distance to rank them
    sorted_routes = sorted(valid_routes.items(), key=lambda item: item[1])
    fastest_route_name, fastest_route_distance = sorted_routes[0]

    print("--- Results ---")
    print("The estimated distances for the routes are:")
    for name, dist in sorted_routes:
        print(f"  - Route {name}: {dist:.0f} meters")

    print(f"\nThe shortest path is Route {fastest_route_name}, making it the fastest option.")
