import math

# Step 1: Define the data provided in the problem.
# The angle of Polaris is equivalent to the latitude.
initial_latitude = 45.0  # Latitude in 1000 CE
latitude_after_first_walk = 44.8 # Latitude in 1100 CE
final_latitude = 34.4 # Latitude in 2000 CE
first_walk_distance_m = 100.0 # Distance of the first walk in meters

# --- Direction Calculation ---
print("--- Step 1: Determine the Direction of Travel ---")
print("The angle of Polaris represents the latitude. The latitude changes from {}° to {}°.".format(initial_latitude, final_latitude))
print("Since the latitude is decreasing, the tree is traveling South.")
print("\nDirection: South\n")


# --- Conversion Factor Calculation ---
print("--- Step 2: Calculate the Distance Conversion Factor ---")
# The first walk is from 1000 CE to 1100 CE.
latitude_change_first_walk = initial_latitude - latitude_after_first_walk
# Convert the first walk's distance to kilometers.
first_walk_distance_km = first_walk_distance_m / 1000.0
# Calculate the conversion factor in km per degree.
km_per_degree = first_walk_distance_km / latitude_change_first_walk
print("The first walk resulted in a latitude change of {:.2f}° ({:.1f}° - {:.1f}°).".format(latitude_change_first_walk, initial_latitude, latitude_after_first_walk))
print("This change corresponds to a known distance of {:.0f} meters, or {:.3f} km.".format(first_walk_distance_m, first_walk_distance_km))
print("Thus, the conversion factor is {:.3f} km / {:.2f}° = {:.1f} km per degree.".format(first_walk_distance_km, latitude_change_first_walk, km_per_degree))


# --- Total Distance Calculation ---
print("\n--- Step 3: Calculate the Total Distance Traveled ---")
# Calculate the total change in latitude over 1000 years.
total_latitude_change = initial_latitude - final_latitude
# Calculate the total distance in kilometers.
total_distance_km = total_latitude_change * km_per_degree
print("The total change in latitude is {:.1f}° ({:.1f}° - {:.1f}°).".format(total_latitude_change, initial_latitude, final_latitude))
print("Using the conversion factor, the total distance is calculated as:")
print("Total Distance (km) = {:.1f} (degrees) * {:.1f} (km/degree) = {:.1f} km".format(total_latitude_change, km_per_degree, total_distance_km))


# --- Final Answer Formatting ---
print("\n--- Step 4: Format the Final Answer ---")
# Apply the final formula as requested.
final_value = total_distance_km * 10
final_answer = round(final_value)
print("The final answer is calculated as: Nearest Integer(Total Distance in km * 10)")
print("Final Calculation: round({:.1f} * 10) = round({}) = {}".format(total_distance_km, int(final_value), final_answer))

print(f"\n<<<{final_answer}>>>")