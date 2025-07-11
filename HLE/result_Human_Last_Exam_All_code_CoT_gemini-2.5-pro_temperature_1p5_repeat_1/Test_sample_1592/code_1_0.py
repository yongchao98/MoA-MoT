import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Earth radius in kilometers
    R = 6371.0

    # Convert coordinates from degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Difference in coordinates
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    # Haversine formula calculation
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c

    return distance

# Coordinates for the Waskaganish gathering place, Quebec
waskaganish_lat = 51.1947
waskaganish_lon = -78.7708

# The Quebec-Ontario border along the James Bay coast is at approximately longitude 79°31′ W.
# We will use a point on this border at the same latitude as Waskaganish for our Ontario reference.
ontario_border_lat = 51.1947
ontario_border_lon = -79.5167

# The closest part of Nunavut is Akimiski Island. We will use a point on its southeastern shore.
nunavut_akimiksi_lat = 52.9
nunavut_akimiksi_lon = -80.9

# Calculate the distances
distance_to_ontario = haversine(waskaganish_lat, waskaganish_lon, ontario_border_lat, ontario_border_lon)
distance_to_nunavut = haversine(waskaganish_lat, waskaganish_lon, nunavut_akimiksi_lat, nunavut_akimiksi_lon)

print("Calculating the closest province or territory to Waskaganish (outside of Quebec).")
print(f"Waskaganish, QC Coordinates: Latitude={waskaganish_lat}, Longitude={waskaganish_lon}")
print("-" * 50)

# Print Ontario calculation details
print("Candidate 1: Ontario")
print(f"Representative Point (ON border): Latitude={ontario_border_lat}, Longitude={ontario_border_lon}")
print(f"Final Equation (Distance in km): haversine({waskaganish_lat}, {waskaganish_lon}, {ontario_border_lat}, {ontario_border_lon})")
print(f"Result: Distance to Ontario is approximately {distance_to_ontario:.2f} km.")
print("-" * 50)

# Print Nunavut calculation details
print("Candidate 2: Nunavut")
print(f"Representative Point (Akimiski Island, NU): Latitude={nunavut_akimiksi_lat}, Longitude={nunavut_akimiksi_lon}")
print(f"Final Equation (Distance in km): haversine({waskaganish_lat}, {waskaganish_lon}, {nunavut_akimiksi_lat}, {nunavut_akimiksi_lon})")
print(f"Result: Distance to Nunavut is approximately {distance_to_nunavut:.2f} km.")
print("-" * 50)


# Determine the closest and print the final answer
if distance_to_ontario < distance_to_nunavut:
    closest_province_territory = "Ontario"
else:
    closest_province_territory = "Nunavut"

print(f"By comparing the two distances, {closest_province_territory} is closer.")
print(f"The closest province or territory to Waskaganish, outside of Quebec, is {closest_province_territory}.")
