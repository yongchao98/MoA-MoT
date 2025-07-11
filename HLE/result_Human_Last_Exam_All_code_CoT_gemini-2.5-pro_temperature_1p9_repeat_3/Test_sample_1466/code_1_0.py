import requests

def get_walking_duration_minutes(waypoints):
    """
    Calculates the walking duration between a series of waypoints using the OSRM API.
    
    Args:
        waypoints (list): A list of (longitude, latitude) tuples.
    
    Returns:
        float: The total walking duration in minutes, or None if the request fails.
    """
    # Format the coordinates into a string like "lon1,lat1;lon2,lat2;..."
    coords_str = ";".join([f"{lon},{lat}" for lon, lat in waypoints])
    
    # Construct the OSRM API URL for the 'foot' profile
    url = f"http://router.project-osrm.org/route/v1/foot/{coords_str}?overview=false"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        data = response.json()
        
        # OSRM returns duration in seconds
        duration_seconds = data['routes'][0]['duration']
        duration_minutes = duration_seconds / 60.0
        return duration_minutes
    except requests.exceptions.RequestException as e:
        print(f"An error occurred during the API request: {e}")
        return None
    except (KeyError, IndexError):
        print("Could not parse the API response.")
        return None

def main():
    """
    Main function to define routes, calculate their durations, and find the fastest one.
    """
    # Define start, end, and intermediate waypoints for each route.
    # Coordinates are in (longitude, latitude) format.
    routes = {
        'A': {
            'description': "via Gresham St and Cheapside (West)",
            'waypoints': [
                (-0.0919, 51.5155),  # Guildhall
                (-0.0954, 51.5150),  # Foster Ln / Cheapside
                (-0.0984, 51.5138)   # St Paul's Cathedral
            ]
        },
        'B': {
            'description': "via London Wall and St Martin's Le Grand",
            'waypoints': [
                (-0.0919, 51.5155),  # Guildhall
                (-0.0953, 51.5179),  # Rotunda / Aldersgate St
                (-0.0970, 51.5155),  # St Martin's Le Grand
                (-0.0984, 51.5138)   # St Paul's Cathedral
            ]
        },
        'C': {
            'description': "via Queen Victoria St",
            'waypoints': [
                (-0.0919, 51.5155),  # Guildhall
                (-0.0886, 51.5134),  # Princes St / Mansion House St
                (-0.0950, 51.5126),  # Cannon St / Queen Victoria St
                (-0.0984, 51.5138)   # St Paul's Cathedral
            ]
        },
        'D': {
            'description': "via Cannon St and Upper Thames St (convoluted)",
            'waypoints': [
                (-0.0919, 51.5155),  # Guildhall
                (-0.0872, 51.5113),  # Cannon St / St Swithin's Ln
                (-0.0977, 51.5106),  # High Timber St
                (-0.0984, 51.5138)   # St Paul's Cathedral
            ]
        },
        'E': {
            'description': "via Barbican and Aldersgate St",
            'waypoints': [
                (-0.0919, 51.5155),  # Guildhall
                (-0.0950, 51.5173),  # Wood St / St Giles Terrace
                (-0.0967, 51.5160),  # Aldersgate St / St Martin's
                (-0.0984, 51.5138)   # St Paul's Cathedral
            ]
        }
    }
    
    results = {}
    print("Calculating walking times for each route...\n")
    
    for route_id, route_data in routes.items():
        duration = get_walking_duration_minutes(route_data['waypoints'])
        if duration is not None:
            results[route_id] = duration
            # Print each number (the calculated time) in the final result
            print(f"Route {route_id} ({route_data['description']}): {duration:.2f} minutes")
        else:
            print(f"Could not calculate time for Route {route_id}.")

    if not results:
        print("\nCould not determine the fastest route.")
        return

    # Find the fastest route
    fastest_route_id = min(results, key=results.get)
    fastest_time = results[fastest_route_id]
    
    print(f"\nBased on the calculations, the fastest route is Route {fastest_route_id} with an estimated time of {fastest_time:.2f} minutes.")
    print(f"\n<<< {fastest_route_id} >>>")


if __name__ == "__main__":
    main()