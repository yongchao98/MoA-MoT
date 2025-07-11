import math

def calculate_distance(lat1, lon1, lat2, lon2):
    """
    Calculates the great-circle distance between two points
    on the earth (specified in decimal degrees).
    """
    # Earth's mean radius in kilometers
    R = 6371.0

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
    Main function to find the closest province/territory to Waskaganish.
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish_lat = 51.1931
    waskaganish_lon = -78.7592

    # A dictionary of the closest points on the borders of neighboring provinces/territories
    # Quebec-Ontario border runs south from James Bay at approx. 79.52° W longitude.
    # The closest part of Nunavut is Akimiski Island in James Bay.
    locations = {
        "Ontario": (51.1931, -79.52),          # Point on the border at the same latitude
        "Nunavut": (52.4, -80.8),              # Southern coast of Akimiski Island
        "Newfoundland and Labrador": (52.2, -67.4), # Point on the western border
        "New Brunswick": (47.8, -69.2)         # Point on the northern border
    }

    distances = {}
    print(f"Finding the closest province or territory to Waskaganish, Quebec (51.19° N, 78.76° W)...")
    print("-" * 30)

    # Calculate and store the distance to each location
    for province, (lat, lon) in locations.items():
        dist = calculate_distance(waskaganish_lat, waskaganish_lon, lat, lon)
        distances[province] = dist
        # Print each number used in the final comparison
        print(f"Distance to the border of {province:<25} is {dist:>7.2f} km.")

    # Find the province with the minimum distance
    closest_province = min(distances, key=distances.get)

    print("-" * 30)
    print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest_province}.")

if __name__ == "__main__":
    main()