import math

# Step 1: Define the given data from the problem description.
initial_latitude_1000CE = 45.0
first_walk_latitude_1100CE = 44.8
final_latitude_2000CE = 34.4

# Distance for the first walk in meters, then converted to kilometers.
first_walk_distance_m = 100.0
first_walk_distance_km = first_walk_distance_m / 1000.0

# --- Calculations ---

# Step 2: Determine the direction of travel.
# In the Northern Hemisphere, a decreasing latitude means moving towards the equator, which is South.
direction = "South"

# Step 3: Calculate the conversion factor (km per degree) using data from the first walk.
latitude_change_first_walk = initial_latitude_1000CE - first_walk_latitude_1100CE
km_per_degree = first_walk_distance_km / latitude_change_first_walk

# Step 4: Calculate the total latitude change and total distance traveled over 1000 years.
total_latitude_change = initial_latitude_1000CE - final_latitude_2000CE
total_distance_km = total_latitude_change * km_per_degree

# --- Output the results ---

print("--- The Ancient Walking Tree's Journey Analysis ---")

print("\n1. Direction of Travel:")
print(f"The tree's latitude (angle to Polaris) decreased from {initial_latitude_1000CE}° to {final_latitude_2000CE}°. This means the tree is consistently walking South.")

print("\n2. Total Distance Traveled:")
print("First, we find the relationship between distance and latitude change from the first walk:")
print(f"  - Latitude change: {initial_latitude_1000CE}° - {first_walk_latitude_1100CE}° = {latitude_change_first_walk:.1f}°")
print(f"  - Distance walked: {first_walk_distance_m} m or {first_walk_distance_km} km")
print(f"This gives us a conversion factor of {km_per_degree:.1f} km per degree of latitude change.")
print("\nNow, we calculate the total journey:")
print(f"  - Total latitude change: {initial_latitude_1000CE}° - {final_latitude_2000CE}° = {total_latitude_change:.1f}°")
print(f"  - Total distance equation: Total Latitude Change * km/degree")
print(f"  - Total distance calculation: {total_latitude_change:.1f}° * {km_per_degree:.1f} km/° = {total_distance_km:.1f} km")

print("\n--- Summary ---")
print(f"Approximate total distance traveled: {total_distance_km:.1f} km")
print(f"Direction of travel: {direction}")

# Step 5: Calculate the final answer in the requested format.
final_answer_value = total_distance_km * 10
final_answer_int = round(final_answer_value)

print("\n--- Final Answer Formulation ---")
print(f"The required value is: Nearest Integer(Total Distance in km * 10)")
print(f"Calculation: round({total_distance_km:.1f} * 10)")
print(f"Calculation: round({final_answer_value})")
print(f"Final Answer: {final_answer_int}")
