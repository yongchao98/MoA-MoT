import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    # Radius of earth in kilometers. Use 6371 for kilometers
    r = 6371
    return c * r

# Coordinates
# Waskaganish, Quebec is located at approximately 51.48° N, 78.76° W
waskaganish_lat, waskaganish_lon = 51.48, -78.76

# The nearest point in Ontario is across James Bay, for example, the community of Moosonee
ontario_point_name = "Moosonee, Ontario"
ontario_lat, ontario_lon = 51.27, -80.64

# The nearest point in Nunavut is Akimiski Island in James Bay
nunavut_point_name = "Akimiski Island, Nunavut"
nunavut_lat, nunavut_lon = 53.00, -81.33


# Calculate distances
dist_to_ontario = haversine(waskaganish_lat, waskaganish_lon, ontario_lat, ontario_lon)
dist_to_nunavut = haversine(waskaganish_lat, waskaganish_lon, nunavut_lat, nunavut_lon)

print(f"Finding the closest province or territory to Waskaganish, QC...")
print(f"Coordinates for Waskaganish, QC: ({waskaganish_lat:.2f} N, {abs(waskaganish_lon):.2f} W)")
print("-" * 30)
print(f"Distance to {ontario_point_name}: {dist_to_ontario:.2f} km")
print(f"Distance to {nunavut_point_name}: {dist_to_nunavut:.2f} km")
print("-" * 30)

if dist_to_ontario < dist_to_nunavut:
    print(f"{ontario_point_name} is closer than {nunavut_point_name}.")
    closest_province_territory = "Ontario"
else:
    print(f"{nunavut_point_name} is closer than {ontario_point_name}.")
    closest_province_territory = "Nunavut"

print(f"\nTherefore, the closest province or territory to Waskaganish, outside of Quebec, is {closest_province_territory}.")
<<<Ontario>>>