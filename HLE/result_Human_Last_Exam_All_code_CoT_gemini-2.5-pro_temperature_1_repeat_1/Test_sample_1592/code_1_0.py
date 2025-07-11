import geopy.distance

def find_closest_province():
    """
    Calculates and prints the closest province or territory to Waskaganish, Quebec.
    """
    # Coordinates for Waskaganish, Quebec, located on Rupert Bay, a part of James Bay.
    waskaganish_coords = (51.483, -78.761)

    # Define coordinates for the closest points of neighboring provinces/territories
    # 1. Ontario: The Quebec-Ontario border south of James Bay is along longitude 79°31′W (-79.517°).
    #    The closest point is on this line at the same latitude as Waskaganish.
    ontario_coords = (51.483, -79.517)

    # 2. Nunavut: The closest part of Nunavut is the southern tip of Akimiski Island in James Bay.
    nunavut_coords = (52.95, -81.33)

    # 3. Newfoundland and Labrador: The westernmost point of its border with Quebec.
    nl_coords = (52.8, -67.1)
    
    # 4. New Brunswick: The northernmost point of its border with Quebec.
    nb_coords = (48.0, -66.9)

    # Calculate the geodesic distance in kilometers for each candidate
    dist_to_ontario = geopy.distance.geodesic(waskaganish_coords, ontario_coords).km
    dist_to_nunavut = geopy.distance.geodesic(waskaganish_coords, nunavut_coords).km
    dist_to_nl = geopy.distance.geodesic(waskaganish_coords, nl_coords).km
    dist_to_nb = geopy.distance.geodesic(waskaganish_coords, nb_coords).km

    distances = {
        "Ontario": dist_to_ontario,
        "Nunavut": dist_to_nunavut,
        "Newfoundland and Labrador": dist_to_nl,
        "New Brunswick": dist_to_nb
    }

    # Find the minimum distance and the corresponding province/territory
    closest_name = min(distances, key=distances.get)
    
    print("This script finds the closest province/territory to Waskaganish, QC.")
    print("-" * 60)
    print(f"Coordinates used for Waskaganish, QC: {waskaganish_coords}")
    print("\nDistances to the nearest border points:")
    print(f"1. To Ontario: {distances['Ontario']:.2f} km")
    print(f"2. To Nunavut: {distances['Nunavut']:.2f} km")
    print(f"3. To Newfoundland and Labrador: {distances['Newfoundland and Labrador']:.2f} km")
    print(f"4. To New Brunswick: {distances['New Brunswick']:.2f} km")
    print("-" * 60)
    print(f"\nThe closest province or territory to Waskaganish, outside of Quebec, is {closest_name}.")

if __name__ == '__main__':
    find_closest_province()