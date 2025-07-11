import geopy.distance

def find_closest_town():
    """
    Finds the closest town with a population over 5,000 to the given house location.
    """
    # Coordinates for the house at 435 N 800 E Rd, Gibson City, IL
    house_coords = (40.4705, -88.2323)

    # Dictionary of nearby towns with population > 5000 and their coordinates
    towns = {
        "Rantoul": (40.3092, -88.1564),
        "Champaign": (40.1164, -88.2434),
        "Bloomington": (40.4842, -88.9937)
    }

    closest_town = None
    min_distance = float('inf')

    print("Calculating distances from the house to nearby towns...")
    for town, coords in towns.items():
        distance = geopy.distance.geodesic(house_coords, coords).miles
        print(f"Distance to {town}: {distance:.2f} miles")
        if distance < min_distance:
            min_distance = distance
            closest_town = town

    print("\nThe closest town with a population over 5,000 is:")
    print(closest_town)

if __name__ == "__main__":
    find_closest_town()