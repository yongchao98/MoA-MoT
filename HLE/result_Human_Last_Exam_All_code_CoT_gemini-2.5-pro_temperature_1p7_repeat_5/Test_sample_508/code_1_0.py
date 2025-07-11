import math

# Step 1 & 2: Define initial parameters and the new location
# Initial sound levels (L1) at 1 meter (r1)
source_levels = {
    'dog': 55,
    'train': 110,
    'construction': 90,
    'people': 75
}

# Source positions
source_positions = {
    'dog': (-25, 0),
    'train': (50, 0),
    'construction': (0, 75),
    'people': (0, -10)
}

# My new position
my_pos = (0, 25)

# Dictionary to hold the calculated sound levels at the new location
levels_at_new_pos = {}

# Step 3 & 4: Calculate new distance and new sound level for each source
print("Calculating sound level from each source at the new location:\n")
for source, pos in source_positions.items():
    # Calculate distance from my new position to the source
    distance = math.sqrt((my_pos[0] - pos[0])**2 + (my_pos[1] - pos[1])**2)
    
    # Get the source's sound level at 1 meter
    l1 = source_levels[source]
    
    # Calculate the sound level (L2) at the new distance
    # Formula: L2 = L1 - 20 * log10(r2/r1), where r1=1
    l2 = l1 - 20 * math.log10(distance)
    
    levels_at_new_pos[source] = l2
    print(f"- {source.capitalize()}:")
    print(f"  - Distance = {distance:.2f} m")
    print(f"  - Sound Level = {l1} - 20*log10({distance:.2f}) = {l2:.2f} dB\n")

# Step 5: Combine the sound levels
# Convert each dB level to intensity, sum them, and convert back to dB
# Formula: L_total = 10 * log10(sum(10^(Li/10)))
sum_of_intensities = 0
for level in levels_at_new_pos.values():
    sum_of_intensities += 10**(level / 10)

total_sound_level = 10 * math.log10(sum_of_intensities)

# Final Output
print("---")
print("To find the total sound level, we use the formula:")
print("L_total = 10 * log10(10^(L1/10) + 10^(L2/10) + ...)\n")
print("The final combined sound level is:\n")

# Retrieve individual levels for the final print statement
l_dog = levels_at_new_pos['dog']
l_train = levels_at_new_pos['train']
l_construction = levels_at_new_pos['construction']
l_people = levels_at_new_pos['people']

# Print the full equation with calculated values
print(f"Total Sound Level = 10 * log10(10^({l_dog:.2f}/10) + 10^({l_train:.2f}/10) + 10^({l_construction:.2f}/10) + 10^({l_people:.2f}/10))")
print(f"= {total_sound_level:.2f} dB")
<<<75.12>>>