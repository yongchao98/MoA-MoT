import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance in kilometers between two points
    on the earth (specified in decimal degrees).
    """
    # Earth radius in kilometers
    R = 6371.0

    # convert decimal degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # differences
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    # haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c

    return distance

def find_closest_province_or_territory():
    """
    Finds the closest province or territory to Waskaganish, QC (outside of Quebec).
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish_coords = {"lat": 51.49, "lon": -78.76}

    # Representative points for neighboring provinces and territories
    locations = {
        "Ontario": {"lat": 51.49, "lon": -79.52},  # Point on the ON-QC border at same latitude
        "Nunavut": {"lat": 52.9, "lon": -81.0},    # Southern coast of Akimiski Island
        "Newfoundland and Labrador": {"lat": 51.42, "lon": -57.11}, # Near the border at Blanc-Sablon
        "New Brunswick": {"lat": 47.37, "lon": -68.32} # City of Edmundston near the border
    }

    distances = {}

    print(f"Calculating distance from Waskaganish, QC ({waskaganish_coords['lat']}° N, {abs(waskaganish_coords['lon'])}° W):\n")

    for loc, coords in locations.items():
        dist = haversine_distance(
            waskaganish_coords["lat"], waskaganish_coords["lon"],
            coords["lat"], coords["lon"]
        )
        distances[loc] = dist
        print(f"To {loc}: approximately {dist:.2f} km")

    # Find the location with the minimum distance
    closest_location = min(distances, key=distances.get)
    min_distance = distances[closest_location]

    print(f"\nThe closest province or territory to Waskaganish, outside of Quebec, is {closest_location}.")

if __name__ == '__main__':
    find_closest_province_or_territory()