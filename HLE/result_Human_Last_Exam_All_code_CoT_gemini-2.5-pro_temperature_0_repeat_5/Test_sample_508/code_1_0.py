import math

def calculate_sound_level_at_distance(L1, d1, d2):
    """Calculates the sound pressure level (L2) at distance d2,
       given the level L1 at distance d1."""
    if d2 <= 0:
        # Cannot be at or inside the sound source
        return -float('inf')
    return L1 - 20 * math.log10(d2 / d1)

def calculate_total_sound_level(levels):
    """Combines multiple sound levels (in dB) to find the total dB."""
    total_intensity = sum(10**(level / 10) for level in levels)
    if total_intensity == 0:
        return -float('inf')
    return 10 * math.log10(total_intensity)

def calculate_distance(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# --- Initial Setup ---
# Sound sources: {name: (L1 @ 1m, (x, y))}
sources = {
    "dog": {"level": 55, "pos": (-25, 0)},
    "train": {"level": 110, "pos": (50, 0)},
    "construction": {"level": 90, "pos": (0, 75)},
    "people": {"level": 75, "pos": (0, -10)}
}

# My new position after walking 25m towards the construction
my_new_pos = (0, 25)
ref_distance = 1.0

# --- Calculations ---
new_levels = {}
for name, data in sources.items():
    distance = calculate_distance(data["pos"], my_new_pos)
    new_level = calculate_sound_level_at_distance(data["level"], ref_distance, distance)
    new_levels[name] = new_level

# Individual sound levels at the new location
L_dog = new_levels["dog"]
L_train = new_levels["train"]
L_construction = new_levels["construction"]
L_people = new_levels["people"]

# Total combined sound level
total_L = calculate_total_sound_level(new_levels.values())

# --- Output ---
print("The sound level from each source at your new location is:")
print(f"Dog: {L_dog:.2f} dB")
print(f"Train: {L_train:.2f} dB")
print(f"Construction: {L_construction:.2f} dB")
print(f"People: {L_people:.2f} dB")
print("\nTo find the total sound level, we combine the intensities:")
print(f"Total dB = 10 * log10(10^({L_dog:.2f}/10) + 10^({L_train:.2f}/10) + 10^({L_construction:.2f}/10) + 10^({L_people:.2f}/10))")
print(f"\nThe total sound level you hear is: {total_L:.2f} dB")
