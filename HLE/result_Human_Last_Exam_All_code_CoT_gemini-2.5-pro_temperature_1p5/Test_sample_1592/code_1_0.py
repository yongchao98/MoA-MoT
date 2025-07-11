import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the distance between two points on Earth using the Haversine formula.
    """
    R = 6371  # Radius of Earth in kilometers
    
    # Convert latitude and longitude from degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    
    # Differences in coordinates
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    
    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    
    return distance

def main():
    """
    Find the closest province/territory to Waskaganish, Quebec.
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish_coords = {"name": "Waskaganish, QC", "lat": 51.49, "lon": -78.76}

    # Coordinates for representative points in neighboring provinces/territories
    locations = {
        "Ontario": {"name": "Moosonee, ON", "lat": 51.27, "lon": -80.64},
        "Nunavut": {"name": "Akimiski Island, NU", "lat": 53.00, "lon": -81.33},
        "Newfoundland and Labrador": {"name": "Labrador City, NL", "lat": 52.94, "lon": -66.91},
        "New Brunswick": {"name": "Edmundston, NB", "lat": 47.37, "lon": -68.32}
    }
    
    # Calculate distances and find the minimum
    closest_province = None
    min_distance = float('inf')

    print(f"Calculating distances from {waskaganish_coords['name']} ({waskaganish_coords['lat']}° N, {waskaganish_coords['lon']}° W):\n")

    for province, coords in locations.items():
        distance = haversine_distance(
            waskaganish_coords['lat'], waskaganish_coords['lon'],
            coords['lat'], coords['lon']
        )
        print(f"Distance to {coords['name']}: {distance:.2f} km")
        
        if distance < min_distance:
            min_distance = distance
            closest_province = province

    print(f"\nThe closest province/territory to Waskaganish, outside of Quebec, is {closest_province}.")

if __name__ == "__main__":
    main()