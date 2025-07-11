import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    # Radius of earth in kilometers. Use 6371 for kilometers
    r = 6371
    return c * r

def find_closest_province():
    """
    Finds the closest province or territory to Waskaganish, QC (excluding Quebec).
    """
    # Coordinates for Waskaganish, Quebec
    waskaganish_lat = 51.48
    waskaganish_lon = -78.76

    # Coordinates for representative cities in other provinces and territories
    # Note: Using Moosonee for Ontario as it is geographically much closer than Toronto.
    locations = {
        "Ontario": (51.27, -80.64), # Moosonee
        "New Brunswick": (45.96, -66.64), # Fredericton
        "Nova Scotia": (44.64, -63.57), # Halifax
        "Manitoba": (49.89, -97.13), # Winnipeg
        "British Columbia": (48.42, -123.36), # Victoria
        "Prince Edward Island": (46.23, -63.12), # Charlottetown
        "Saskatchewan": (50.44, -104.61), # Regina
        "Alberta": (53.54, -113.49), # Edmonton
        "Newfoundland and Labrador": (47.56, -52.71), # St. John's
        "Northwest Territories": (62.45, -114.37), # Yellowknife
        "Yukon": (60.72, -135.05), # Whitehorse
        "Nunavut": (63.75, -68.52) # Iqaluit
    }

    closest_province = None
    min_distance = float('inf')

    for province, (lat, lon) in locations.items():
        distance = haversine(waskaganish_lat, waskaganish_lon, lat, lon)
        if distance < min_distance:
            min_distance = distance
            closest_province = province

    print(f"The Waskaganish gathering place is located at latitude {waskaganish_lat} and longitude {waskaganish_lon}.")
    print("Calculating distances to other provinces and territories...")
    print(f"The closest province/territory to Waskaganish (outside of Quebec) is {closest_province}.")
    print(f"The approximate distance to the nearest point in {closest_province} (Moosonee) is {min_distance:.2f} km.")

if __name__ == '__main__':
    find_closest_province()