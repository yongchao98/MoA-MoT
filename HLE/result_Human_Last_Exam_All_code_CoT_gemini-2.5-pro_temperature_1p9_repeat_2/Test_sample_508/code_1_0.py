import math

def calculate_new_sound_level(L1, r2):
    """
    Calculates the sound level (L2) at a new distance (r2),
    given the sound level (L1) at 1 meter.
    Formula: L2 = L1 - 20 * log10(r2)
    """
    if r2 <= 0:
        return float('inf') # Cannot be at or inside the source
    return L1 - 20 * math.log10(r2)

def calculate_distance(p1, p2):
    """Calculates the Euclidean distance between two points (x, y)."""
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Step 1 & 2: Define initial positions and new position
my_new_pos = (0, 25)

sources = {
    "Dog": {"pos": (-25, 0), "L1": 55},
    "Train": {"pos": (50, 0), "L1": 110},
    "Construction": {"pos": (0, 75), "L1": 90},
    "People": {"pos": (0, -10), "L1": 75}
}

print("Step 1: Calculating the new sound level from each source at your new position.\n")

# Step 3 & 4: Calculate new distances and sound levels for each source
new_levels_dB = {}
intensity_terms = {}

for name, data in sources.items():
    distance = calculate_distance(data["pos"], my_new_pos)
    new_dB = calculate_new_sound_level(data["L1"], distance)
    new_levels_dB[name] = new_dB
    intensity_terms[name] = 10**(new_dB / 10)
    print(f"Source: {name}")
    print(f"  - Distance from new position: {distance:.2f} meters")
    print(f"  - Sound Level = {data['L1']} - 20 * log10({distance:.2f}) = {new_dB:.2f} dB\n")

# Step 5: Calculate total sound level
total_intensity_term = sum(intensity_terms.values())
total_dB = 10 * math.log10(total_intensity_term)

print("Step 2: Calculating the total sound level by summing the relative intensities.\n")
print("Total Sound Level Equation:")
print(f"L_total = 10 * log10( 10^({new_levels_dB['Dog']:.2f}/10) + 10^({new_levels_dB['Train']:.2f}/10) + 10^({new_levels_dB['Construction']:.2f}/10) + 10^({new_levels_dB['People']:.2f}/10) )")
print("\nPlugging in the intensity values:")
print(f"L_total = 10 * log10( {intensity_terms['Dog']:.2f} + {intensity_terms['Train']:.2f} + {intensity_terms['Construction']:.2f} + {intensity_terms['People']:.2f} )")
print(f"L_total = 10 * log10( {total_intensity_term:.2f} )")
print("\nFinal Result:")
print(f"The total sound level you hear is {total_dB:.2f} dB.")
print(f"<<<{total_dB:.2f}>>>")