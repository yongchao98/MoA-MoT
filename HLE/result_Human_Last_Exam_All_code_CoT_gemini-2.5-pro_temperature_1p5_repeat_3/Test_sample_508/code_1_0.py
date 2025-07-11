import math

# --- 1. Define initial parameters ---

# Sound sources: {name: [level_at_1m_dB, x_pos, y_pos]}
sources = {
    'dog': [55, -25, 0],
    'train': [110, 50, 0],
    'construction': [90, 0, 75],
    'people': [75, 0, -10]
}

# Your new position after walking 25m towards the construction
# Initial position was (0,0), walking towards (0, 75)
my_new_pos = (0, 25)

# --- 2. Calculate sound levels at the new position ---

sound_levels_at_new_pos_dB = {}
sound_intensities = []

print("Calculating the sound level from each source at your new location:\n")

for name, params in sources.items():
    level_at_1m = params[0]
    source_pos_x = params[1]
    source_pos_y = params[2]
    
    # Calculate new distance using the distance formula
    # d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    new_distance = math.sqrt((source_pos_x - my_new_pos[0])**2 + (source_pos_y - my_new_pos[1])**2)
    
    # Calculate the sound level (dB) at the new distance
    # L2 = L1 - 20 * log10(d2 / d1), where d1=1m
    level_at_new_distance = level_at_1m - 20 * math.log10(new_distance)
    sound_levels_at_new_pos_dB[name] = level_at_new_distance
    
    # Convert dB to a relative intensity for summing
    # I = 10^(L/10)
    intensity = 10**(level_at_new_distance / 10)
    sound_intensities.append(intensity)
    
    print(f"Source: {name.capitalize()}")
    print(f"  - New distance: {new_distance:.2f} meters")
    print(f"  - Sound level at new distance: {level_at_new_distance:.2f} dB\n")

# --- 3. Sum intensities and convert back to total dB ---

# Sum the individual intensities
total_intensity = sum(sound_intensities)

# Convert total intensity back to decibels
# L_total = 10 * log10(I_total)
total_sound_level_dB = 10 * math.log10(total_intensity)

# --- 4. Print the final result and equation ---

# Retrieve the calculated dB levels for the final equation display
L_dog = sound_levels_at_new_pos_dB['dog']
L_train = sound_levels_at_new_pos_dB['train']
L_construction = sound_levels_at_new_pos_dB['construction']
L_people = sound_levels_at_new_pos_dB['people']

print("--- Final Calculation ---")
print("To find the total sound level, we sum the sound intensities and convert back to dB.")
print("\nFinal Equation:")
print(f"Total Sound Level = 10 * log10(10^({L_dog:.2f}/10) + 10^({L_train:.2f}/10) + 10^({L_construction:.2f}/10) + 10^({L_people:.2f}/10))")

print("\nResult:")
print(f"The total sound level you hear is: {total_sound_level_dB:.2f} dB")
<<<75.12>>>