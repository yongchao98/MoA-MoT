import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Earth's radius in kilometers
    R = 6371.0
    
    # Convert decimal degrees to radians
    lat1_rad, lon1_rad, lat2_rad, lon2_rad = map(math.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    distance = R * c
    return distance

def solve():
    """
    Main function to find the closest province/territory to Waskaganish.
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish = {"name": "Waskaganish", "lat": 51.1939, "lon": -78.7614}

    # Coordinates for the closest points in neighboring provinces/territories
    locations = [
        {
            "name": "Ontario",
            "type": "border",
            "lat": 51.463, # Point on the Ontario-Quebec land border meeting James Bay
            "lon": -79.518
        },
        {
            "name": "Nunavut",
            "type": "landmass (Charlton Island)", # Closest significant landmass
            "lat": 52.00, 
            "lon": -79.40
        },
        {
            "name": "Newfoundland and Labrador",
            "type": "border", # Point on the Quebec-Labrador border
            "lat": 52.00,
            "lon": -65.0
        }
    ]

    distances = {}
    print(f"Finding the closest province or territory to Waskaganish, QC ({waskaganish['lat']}° N, {abs(waskaganish['lon'])}° W).\n")

    for loc in locations:
        dist = haversine_distance(waskaganish["lat"], waskaganish["lon"], loc["lat"], loc["lon"])
        distances[loc["name"]] = dist
        
        print(f"Calculating distance to {loc['name']}:")
        
        # Format the numbers for the equation representation
        lat1_str = f"{waskaganish['lat']:.4f}"
        lon1_str = f"{waskaganish['lon']:.4f}"
        lat2_str = f"{loc['lat']:.4f}"
        lon2_str = f"{loc['lon']:.4f}"
        dist_str = f"{dist:.2f}"
        
        # Print the simplified equation with the inputs and output
        print(f"  Start: Waskaganish ({lat1_str}, {lon1_str})")
        print(f"  End:   {loc['name']} ({lat2_str}, {lon2_str})")
        print(f"  haversine({lat1_str}, {lon1_str}, {lat2_str}, {lon2_str}) = {dist_str} km\n")

    # Find the location with the minimum distance
    closest_location_name = min(distances, key=distances.get)
    min_distance = distances[closest_location_name]

    print(f"Comparing the distances, the closest province or territory to Waskaganish, outside of Quebec, is {closest_location_name}.")

solve()
<<<Ontario>>>