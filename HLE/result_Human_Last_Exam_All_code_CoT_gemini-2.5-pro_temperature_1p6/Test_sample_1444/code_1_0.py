import math

# --- Data from the problem ---
# Year: 1000, 1100, ..., 2000
angles = [45.0, 44.8, 44.3, 43.5, 42.4, 41.0, 39.4, 37.8, 36.5, 35.2, 34.4]
initial_angle = angles[0]
angle_after_first_walk = angles[1]
final_angle = angles[-1]
distance_first_walk_m = 100.0

# --- Step 1: Determine the direction ---
# The angle of Polaris approximates latitude. A decrease in latitude in the
# Northern Hemisphere means moving South.
direction = "South"
print(f"Direction of Travel: The tree is walking {direction}.\n")

# --- Step 2: Calibrate the distance scale using the first walk ---
print("--- Calculating the Total Distance ---")
# The change in angle for the first 100-meter walk
angular_change_first_walk = initial_angle - angle_after_first_walk

# The calibration factor (meters per degree)
meters_per_degree = distance_first_walk_m / angular_change_first_walk

print(f"The first walk from 1000 to 1100 CE gives us a scale.")
print(f"An angular change of {initial_angle}° - {angle_after_first_walk}° = {angular_change_first_walk:.1f}° corresponded to {distance_first_walk_m:.0f} meters.")
print(f"This establishes a scale factor of {meters_per_degree:.0f} meters per degree.\n")

# --- Step 3: Calculate the total distance ---
# The total change in angle over the 1000 years
total_angular_change = initial_angle - final_angle

# Calculate the total distance using the scale factor
total_distance_m = total_angular_change * meters_per_degree
total_distance_km = total_distance_m / 1000

print(f"The total angular change from 1000 CE to 2000 CE is {initial_angle}° - {final_angle}° = {total_angular_change:.1f}°.")
print("\nThe final equation for the total distance in meters is:")
print(f"Total Distance (m) = (Total Angular Change) * (Meters Per Degree)")
print(f"Total Distance (m) = {total_angular_change:.1f} * ({distance_first_walk_m:.0f} / {angular_change_first_walk:.1f}) = {total_distance_m:.0f} meters.")

print(f"\nThe total distance traveled over 1000 years is {total_distance_m:.0f} meters, which is {total_distance_km:.1f} kilometers.\n")

# --- Step 4: Format the final answer ---
print("--- Calculating the Final Answer ---")
# Calculate the value as per the requested format
value_for_rounding = total_distance_km * 10
final_answer = round(value_for_rounding)

print(f"The required format is: Nearest Integer(Total Distance in km * 10)")
print(f"Calculation: round({total_distance_km:.1f} * 10) = round({value_for_rounding:.1f}) = {final_answer}")

print(f"<<<{final_answer}>>>")