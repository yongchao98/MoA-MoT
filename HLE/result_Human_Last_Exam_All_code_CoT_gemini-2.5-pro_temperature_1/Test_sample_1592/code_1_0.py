import geopy.distance

def find_closest_province():
    """
    Calculates the closest province or territory to Waskaganish, QC, excluding Quebec itself.
    """
    # Coordinates for Waskaganish, QC
    waskaganish_coords = (51.1950, -78.7608)

    # Coordinates for representative points in neighboring provinces/territories
    locations = {
        "Ontario": (51.1950, -80.5), # A point on the ON coast directly west of Waskaganish
        "Nunavut": (52.7, -79.7), # A point on the southeastern coast of Akimiski Island, NU
        "Newfoundland and Labrador": (52.9463, -66.9133), # Labrador City, NL
    }

    print(f"Finding the closest province/territory to Waskaganish, QC {waskaganish_coords}...\n")

    distances = {}
    for province, coords in locations.items():
        # Calculate the distance using the geodesic formula (most accurate)
        distance_km = geopy.distance.geodesic(waskaganish_coords, coords).km
        distances[province] = distance_km
        print(f"Distance from Waskaganish to {province} at {coords}: {distance_km:.1f} km")

    # Find the province with the minimum distance
    closest_province = min(distances, key=distances.get)
    min_distance = distances[closest_province]

    print(f"\nBased on the calculations, the closest province or territory to Waskaganish is {closest_province}.")

if __name__ == '__main__':
    find_closest_province()
