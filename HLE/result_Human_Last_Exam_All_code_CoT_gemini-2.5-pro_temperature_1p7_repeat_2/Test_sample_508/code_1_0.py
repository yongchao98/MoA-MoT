import math

# Step 1: Define initial conditions for sound sources
# Format: {name: {'pos': (x, y), 'spl1': sound_level_at_1m}}
sources = {
    'dog': {'pos': (-25, 0), 'spl1': 55},
    'train': {'pos': (50, 0), 'spl1': 110},
    'construction': {'pos': (0, 75), 'spl1': 90},
    'people': {'pos': (0, -10), 'spl1': 75}
}

# Step 2: Define your new position after walking 25m towards the construction
my_new_pos = (0, 25)

# Step 3 & 4: Calculate new distances and new Sound Pressure Levels (SPL)
new_spl_values = []
sum_of_intensities = 0

print("Calculating individual sound levels at the new location:\n")

for name, data in sources.items():
    source_pos = data['pos']
    spl1 = data['spl1']
    
    # Step 3: Calculate new distance (r2)
    distance = math.sqrt((my_new_pos[0] - source_pos[0])**2 + (my_new_pos[1] - source_pos[1])**2)
    
    # Step 4: Calculate new SPL (Lp2)
    # Lp2 = Lp1 - 20 * log10(r2/r1), where r1=1m
    spl2 = spl1 - 20 * math.log10(distance)
    new_spl_values.append(spl2)
    
    print(f"Sound level from {name.capitalize()}: {spl2:.1f} dB (at {distance:.1f} meters)")
    
    # Prepare for Step 5: Convert SPL to intensity and add to a running sum
    sum_of_intensities += 10**(spl2 / 10)

# Step 5: Calculate the total combined SPL
total_spl = 10 * math.log10(sum_of_intensities)

# Step 6: Print the final calculation and result
print("\n---")
print("The final calculation combines these decibel levels:")
# Build the equation string using the calculated individual sound levels
equation = (
    f"Total dB = 10 * log10( 10^({new_spl_values[0]:.1f}/10) "
    f"+ 10^({new_spl_values[1]:.1f}/10) "
    f"+ 10^({new_spl_values[2]:.1f}/10) "
    f"+ 10^({new_spl_values[3]:.1f}/10) )"
)
print(equation)
print(f"\nTotal Sound Level = {total_spl:.1f} dB")
<<<75.1>>>