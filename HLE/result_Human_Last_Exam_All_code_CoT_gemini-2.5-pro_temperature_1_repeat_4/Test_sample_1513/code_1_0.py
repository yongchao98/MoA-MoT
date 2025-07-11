import math

def find_island_town():
    """
    This script identifies the island town from the riddle and verifies its remoteness.
    The town is Hanga Roa on Easter Island. We will calculate its distance to the
    nearest inhabited neighbor, Pitcairn Island, to confirm it's over 500 km away.
    """

    # The town name is the solution to the riddle.
    town_name = "Hanga Roa"

    # The distance from the riddle's condition.
    riddle_distance_km = 500

    # Coordinates in decimal degrees (Latitude, Longitude)
    # Hanga Roa, Easter Island
    lat1, lon1 = -27.1436, -109.4267
    # Adamstown, Pitcairn Island
    lat2, lon2 = -25.0667, -130.1000

    # Earth's radius in kilometers
    R = 6371

    # Convert coordinates from degrees to radians
    lat1_rad, lon1_rad = math.radians(lat1), math.radians(lon1)
    lat2_rad, lon2_rad = math.radians(lat2), math.radians(lon2)

    # Difference in coordinates
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    # Haversine formula to calculate the distance
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    calculated_distance_km = R * c

    print(f"The island town described is: {town_name}")
    print(f"To verify its remoteness, we calculate the distance to its nearest inhabited neighbor (Pitcairn Island).")
    print(f"The calculated distance is approximately {int(calculated_distance_km)} kilometers.")
    print(f"This is significantly greater than the {riddle_distance_km} kilometers mentioned in the riddle.")

# Execute the function to print the results
find_island_town()