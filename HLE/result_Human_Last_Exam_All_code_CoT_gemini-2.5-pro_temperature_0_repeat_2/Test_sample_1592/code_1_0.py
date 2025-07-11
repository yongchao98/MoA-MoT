import geopy.distance

def find_closest_province():
    """
    Calculates and determines the closest province or territory to Waskaganish, QC (outside of Quebec).
    """
    # Step 1: Define coordinates for Waskaganish
    waskaganish = (51.473, -78.761)

    # Step 2: Define coordinates for the closest points in neighboring regions
    # The Ontario border is at longitude 79° 31' W (-79.517 degrees).
    # The closest point is on this line at the same latitude as Waskaganish.
    ontario_border_point = (51.473, -79.517)

    # The closest point in Nunavut is the southern tip of Charlton Island in James Bay.
    nunavut_charlton_island = (51.99, -79.42)

    # Step 3: Calculate the distances
    # The geopy.distance.geodesic function calculates the distance between two points on the Earth's surface.
    # We will treat the coordinates as the numbers in our "equation".
    
    print("Equation 1: Distance to Ontario")
    print(f"Calculating geodesic distance between Waskaganish {waskaganish} and the Ontario border {ontario_border_point}")
    dist_to_ontario = geopy.distance.geodesic(waskaganish, ontario_border_point).km
    print(f"Result: {dist_to_ontario:.2f} km\n")

    print("Equation 2: Distance to Nunavut")
    print(f"Calculating geodesic distance between Waskaganish {waskaganish} and Nunavut (Charlton Island) {nunavut_charlton_island}")
    dist_to_nunavut = geopy.distance.geodesic(waskaganish, nunavut_charlton_island).km
    print(f"Result: {dist_to_nunavut:.2f} km\n")

    # Step 4: Compare distances and determine the closest
    print("--- Comparison ---")
    if dist_to_ontario < dist_to_nunavut:
        closest_province = "Ontario"
        closest_distance = dist_to_ontario
        print(f"{dist_to_ontario:.2f} km (to Ontario) < {dist_to_nunavut:.2f} km (to Nunavut)")
    else:
        closest_province = "Nunavut"
        closest_distance = dist_to_nunavut
        print(f"{dist_to_nunavut:.2f} km (to Nunavut) < {dist_to_ontario:.2f} km (to Ontario)")

    print(f"\nThe closest province or territory to Waskaganish, outside of Quebec, is {closest_province}.")

if __name__ == '__main__':
    find_closest_province()