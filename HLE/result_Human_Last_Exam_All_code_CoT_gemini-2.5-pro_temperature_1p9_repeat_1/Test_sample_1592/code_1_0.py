import math

def calculate_and_print_distance(lat1_deg, lon1_deg, lat2_deg, lon2_deg, origin_name, destination_name):
    """
    Calculates the Haversine distance between two points and prints the steps.
    """
    print(f"--- Calculating distance from {origin_name} to {destination_name} ---")

    # Convert decimal degrees to radians
    lon1_rad, lat1_rad, lon2_rad, lat2_rad = map(math.radians, [lon1_deg, lat1_deg, lon2_deg, lat2_deg])
    
    # Earth's radius in kilometers
    R = 6371

    # Haversine formula steps
    # Step 1: Calculate changes in latitude and longitude
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad
    
    # Step 2: Calculate 'a'
    # a = sin²(Δlat/2) + cos(lat1) * cos(lat2) * sin²(Δlon/2)
    a = math.sin(delta_lat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2)**2
    print(f"Equation for 'a': a = sin²({delta_lat:.4f}/2) + cos({lat1_rad:.4f}) * cos({lat2_rad:.4f}) * sin²({delta_lon:.4f}/2)")
    print(f"Result for 'a': a = {a:.6f}")
    
    # Step 3: Calculate 'c'
    # c = 2 * atan2(√a, √(1-a))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    print(f"Equation for 'c': c = 2 * atan2(√{a:.6f}, √(1-{a:.6f}))")
    print(f"Result for 'c': c = {c:.6f}")

    # Step 4: Calculate final distance
    # d = R * c
    distance = R * c
    print(f"Equation for distance: d = {R} * {c:.6f}")
    print(f"Final distance: {distance:.2f} km\n")

    return distance

# Coordinates of the locations
waskaganish_loc = {"name": "Waskaganish, QC", "coords": (51.49, -78.76)}
moosonee_loc = {"name": "Moosonee, ON", "coords": (51.27, -80.64)}
sanikiluaq_loc = {"name": "Sanikiluaq, NU", "coords": (56.54, -79.23)}

# Calculate distance to Ontario candidate
dist_to_ontario = calculate_and_print_distance(
    waskaganish_loc["coords"][0], waskaganish_loc["coords"][1],
    moosonee_loc["coords"][0], moosonee_loc["coords"][1],
    waskaganish_loc["name"], moosonee_loc["name"]
)

# Calculate distance to Nunavut candidate
dist_to_nunavut = calculate_and_print_distance(
    waskaganish_loc["coords"][0], waskaganish_loc["coords"][1],
    sanikiluaq_loc["coords"][0], sanikiluaq_loc["coords"][1],
    waskaganish_loc["name"], sanikiluaq_loc["name"]
)

# Compare distances and print the conclusion
print("--- Conclusion ---")
if dist_to_ontario < dist_to_nunavut:
    closest_province = "Ontario"
    print(f"The distance to Moosonee, Ontario ({dist_to_ontario:.2f} km) is less than the distance to Sanikiluaq, Nunavut ({dist_to_nunavut:.2f} km).")
else:
    closest_province = "Nunavut"
    print(f"The distance to Sanikiluaq, Nunavut ({dist_to_nunavut:.2f} km) is less than the distance to Moosonee, Ontario ({dist_to_ontario:.2f} km).")

print(f"\nThe closest province or territory to Waskaganish, outside of Quebec, is {closest_province}.")

<<<Ontario>>>