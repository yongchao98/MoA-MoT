import math

# Step 1: Define the given positions and the distance of the first walk.
pos_1000 = 45.0  # Initial position in 1000 CE
pos_1100 = 44.8  # Position after first walk
pos_2000 = 34.4  # Final position in 2000 CE
dist_first_walk_m = 100.0  # Distance of the first walk in meters

print("--- Analysis of the Ancient Walking Tree's Journey ---")
print("\n1. Determining the Direction of Travel:")
if pos_1000 > pos_2000:
    direction = "South"
    print(f"The latitude decreased from {pos_1000}° to {pos_2000}°. This indicates the tree is moving South.")
else:
    # This case is not expected based on the data
    direction = "North"
    print(f"The latitude increased from {pos_1000}° to {pos_2000}°. This indicates the tree is moving North.")

print("\n2. Calculating the Distance-to-Latitude Conversion Factor:")
# The first walk's details are used for calibration
lat_change_first_walk = pos_1000 - pos_1100
dist_first_walk_km = dist_first_walk_m / 1000.0

# This factor tells us how many km the tree travels for each degree of latitude change
km_per_degree = dist_first_walk_km / lat_change_first_walk
print(f"The first walk covered {dist_first_walk_m} meters ({dist_first_walk_km} km).")
print(f"The latitude change during this walk was {pos_1000}° - {pos_1100}° = {lat_change_first_walk:.1f}°.")
print(f"This establishes a conversion factor of {km_per_degree:.1f} km per degree of latitude change.")

print("\n3. Calculating the Total Distance Traveled:")
# Total change in latitude over 1000 years
total_lat_change = pos_1000 - pos_2000

# Calculate total distance using the conversion factor
total_distance_km = total_lat_change * km_per_degree

print(f"The total latitude change over 1000 years is {pos_1000}° - {pos_2000}° = {total_lat_change:.1f}°.")
print("Using the conversion factor, we can find the total distance:")
# Print the final equation as requested
print(f"Total Distance = ({pos_1000} - {pos_2000}) * ({dist_first_walk_m} / ({pos_1000} - {pos_1100}) / 1000.0) = {total_distance_km:.1f} km")


print("\n--- Summary ---")
print(f"Direction of Travel: {direction}")
print(f"Approximate Total Distance Traveled: {total_distance_km:.1f} km")

# Calculate and format the final answer as per the user's request
final_answer_value = round(total_distance_km * 10)
print(f"\nFinal answer calculation: Nearest Integer({total_distance_km:.1f} * 10) = Nearest Integer({total_distance_km*10:.1f}) = {final_answer_value}")

print("<<<" + str(final_answer_value) + ">>>")