import math

def calculate_distance(lat1, lon1, lat2, lon2):
    """
    Calculates the great-circle distance between two points
    on the earth (specified in decimal degrees) using the Haversine formula.
    """
    # Earth radius in kilometers
    R = 6371

    # Convert decimal degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Haversine formula
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance

def main():
    """
    Main function to find the closest province/territory to Waskaganish.
    """
    # Coordinates for the Waskaganish gathering place in Quebec
    waskaganish_lat, waskaganish_lon = 51.1944, -78.7583

    # Coordinates for the closest points in neighboring provinces/territories
    # New Brunswick and Newfoundland/Labrador are clearly much farther and are not included in this direct comparison.
    locations = {
        "Ontario": {
            "coords": (51.1944, -79.52),  # A point on the ON/QC border directly west of Waskaganish
            "description": "the Ontario border"
        },
        "Nunavut": {
            "coords": (52.6, -80.8),  # The southeastern coast of Akimiski Island, which is part of Nunavut
            "description": "Akimiski Island (Nunavut)"
        }
    }

    closest_location_name = None
    min_distance = float('inf')

    print(f"Calculating distances from Waskaganish, QC ({waskaganish_lat}, {waskaganish_lon})...\n")

    for name, data in locations.items():
        dest_lat, dest_lon = data["coords"]
        distance = calculate_distance(waskaganish_lat, waskaganish_lon, dest_lat, dest_lon)

        print(f"The distance from Waskaganish to {data['description']} is {distance:.2f} km.")

        if distance < min_distance:
            min_distance = distance
            closest_location_name = name

    print(f"\nBased on the calculations, the closest province or territory to the Waskaganish gathering place, outside of Quebec, is {closest_location_name}.")


if __name__ == "__main__":
    main()