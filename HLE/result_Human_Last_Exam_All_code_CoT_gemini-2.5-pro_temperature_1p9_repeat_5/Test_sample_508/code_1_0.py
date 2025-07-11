import math

# Step 1 & 2: Define initial positions, sound levels, and the new observation point.
sources = {
    "dog": {"pos": (-25, 0), "level_at_1m": 55},
    "train": {"pos": (50, 0), "level_at_1m": 110},
    "construction": {"pos": (0, 75), "level_at_1m": 90},
    "people": {"pos": (0, -10), "level_at_1m": 75}
}
observer_pos = (0, 25)

print("Calculating the total sound level at the new position (0, 25).\n")

# Store the individual levels at the new location
levels_at_new_pos = {}
total_intensity_sum = 0

# Step 3 & 4: Calculate new distance and new sound level for each source
for name, data in sources.items():
    source_pos = data["pos"]
    level_at_1m = data["level_at_1m"]
    
    # Calculate the new distance
    distance = math.sqrt((source_pos[0] - observer_pos[0])**2 + (source_pos[1] - observer_pos[1])**2)
    
    # Calculate the sound level at the new distance
    level_at_distance = level_at_1m - 20 * math.log10(distance)
    levels_at_new_pos[name] = level_at_distance
    
    print(f"The sound level from the {name} (distance: {distance:.2f} m) is {level_at_distance:.2f} dB.")
    
    # Calculate the intensity term for the final summation
    intensity_term = 10**(level_at_distance / 10)
    total_intensity_sum += intensity_term

print("\nTo find the total sound level, we use the formula:")
print("L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + ...)\n")

# Step 5: Combine the sound levels
total_level = 10 * math.log10(total_intensity_sum)

# Print the final detailed equation
level_dog = levels_at_new_pos['dog']
level_train = levels_at_new_pos['train']
level_construction = levels_at_new_pos['construction']
level_people = levels_at_new_pos['people']

print("Final Equation:")
print(f"L_total = 10 * log10( 10^({level_dog:.2f}/10) + 10^({level_train:.2f}/10) + 10^({level_construction:.2f}/10) + 10^({level_people:.2f}/10) )")
print(f"\nThe total sound level you hear is: {total_level:.2f} dB")

<<<85.44>>>