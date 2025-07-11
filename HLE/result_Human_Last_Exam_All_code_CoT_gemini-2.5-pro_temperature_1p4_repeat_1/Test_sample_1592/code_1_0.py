import math

def calculate_and_print_distance(name1, lat1, lon1, name2, lat2, lon2):
    """
    Calculates the Haversine distance between two points and prints the steps.
    """
    R = 6371  # Earth radius in kilometers

    print(f"\n--- Calculating distance between {name1} and the closest point in {name2} ---")
    print(f"Coordinates for {name2}: (lat={lat2}, lon={lon2})")

    # Convert degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Differences in coordinates
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad
    
    # Haversine formula part 1: a = sin²(Δφ/2) + cos(φ1) * cos(φ2) * sin²(Δλ/2)
    a = math.sin(delta_lat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2)**2
    print(f"Equation for 'a': sin²({delta_lat:.4f}/2) + cos({lat1_rad:.4f}) * cos({lat2_rad:.4f}) * sin²({delta_lon:.4f}/2)")
    print(f"Result 'a': {a:.8f}")

    # Haversine formula part 2: c = 2 * atan2(√a, √(1-a))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    print(f"Equation for 'c': 2 * atan2(√{a:.8f}, √(1-{a:.8f}))")
    print(f"Result 'c': {c:.8f}")

    # Haversine formula part 3: d = R * c
    distance = R * c
    print(f"Equation for distance 'd': {R} * {c:.8f}")
    print(f"Final Distance: {distance:.2f} km")

    return distance

def find_closest_province_territory():
    """
    Finds the closest province or territory to Waskaganish, outside of Quebec.
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish_name = "Waskaganish"
    waskaganish_lat = 51.1939
    waskaganish_lon = -78.7606
    print(f"Base location: {waskaganish_name}, Quebec (lat={waskaganish_lat}, lon={waskaganish_lon})")

    # Coordinates for the closest points in neighboring provinces/territories
    locations = {
        "Ontario": (51.1939, -79.5167),  # A point on the border due west of Waskaganish
        "Nunavut": (52.6, -81.5),  # Southern coast of Akimiski Island
        "Newfoundland and Labrador": (52.94, -66.91),  # Labrador City, near the border
        "New Brunswick": (48.05, -67.75)  # A point on the border near Edmundston
    }

    closest_location_name = None
    min_distance = float('inf')

    # Iterate through locations, calculate distance, and find the minimum
    for name, (lat, lon) in locations.items():
        dist = calculate_and_print_distance(waskaganish_name, waskaganish_lat, waskaganish_lon, name, lat, lon)
        if dist < min_distance:
            min_distance = dist
            closest_location_name = name
    
    print("\n--- Conclusion ---")
    print(f"Comparing all distances, the closest province or territory is {closest_location_name} with a distance of {min_distance:.2f} km.")
    
    # Final answer in the specified format
    print(f"\n<<<{closest_location_name}>>>")


if __name__ == '__main__':
    find_closest_province_territory()