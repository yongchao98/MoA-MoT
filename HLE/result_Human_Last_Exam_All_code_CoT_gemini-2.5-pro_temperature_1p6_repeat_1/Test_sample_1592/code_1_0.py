import geopy.distance

def find_closest_province():
    """
    Calculates the closest province or territory to Waskaganish, QC, from a list of candidates.
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish = (51.4886, -78.7694)

    # Coordinates for the nearest points of surrounding provinces/territories
    locations = {
        "Ontario": (51.4886, -79.5167),  # A point on the ON-QC land border due west of Waskaganish
        "Nunavut": (52.88, -80.88),        # Southeast tip of Akimiski Island, Nunavut, in James Bay
        "Newfoundland and Labrador": (52.94, -66.82), # Point on the QC-Labrador border near Labrador City
        "New Brunswick": (48.05, -67.50)     # Northernmost point of the NB-QC border
    }

    distances = {}
    print("Calculating distances from Waskaganish, QC:")
    print("-" * 40)

    # Calculate distance to each location
    for province, coords in locations.items():
        dist = geopy.distance.geodesic(waskaganish, coords).km
        distances[province] = dist
        print(f"Distance to the {province} border: {dist:.2f} km")

    # Find the minimum distance
    closest_province = min(distances, key=distances.get)
    min_distance = distances[closest_province]

    print("-" * 40)
    print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest_province}.")
    print(f"The distance to its border is approximately {min_distance:.2f} km.")

if __name__ == '__main__':
    find_closest_province()