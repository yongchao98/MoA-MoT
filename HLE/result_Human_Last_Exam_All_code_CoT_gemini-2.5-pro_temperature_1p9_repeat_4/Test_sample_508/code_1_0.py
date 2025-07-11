import math

def calculate_new_sound_level(Lp1, r2):
    """
    Calculates the sound pressure level (Lp2) at a new distance (r2)
    from a source, given the level (Lp1) at 1 meter.
    Formula: Lp2 = Lp1 - 20 * log10(r2)
    """
    if r2 <= 0:
        return float('inf') # Cannot be at or inside the source
    return Lp1 - 20 * math.log10(r2)

def calculate_total_sound_level(levels):
    """
    Combines multiple sound levels into a total perceived sound level.
    Formula: L_total = 10 * log10( Σ(10^(Li / 10)) )
    """
    total_intensity = sum(10**(level / 10) for level in levels)
    return 10 * math.log10(total_intensity)

def calculate_distance(pos1, pos2):
    """Calculates the Euclidean distance between two points."""
    return math.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)

# Initial setup
my_initial_pos = (0, 0)
my_new_pos = (0, 25) # Walked 25m towards construction at (0, 75)

sources = {
    "dog": {"Lp1": 55, "pos": (-25, 0)},
    "train": {"Lp1": 110, "pos": (50, 0)},
    "construction": {"Lp1": 90, "pos": (0, 75)},
    "people": {"Lp1": 75, "pos": (0, -10)}
}

# --- Calculations ---

# 1. Calculate new distances to each source
new_distances = {name: calculate_distance(my_new_pos, data["pos"]) for name, data in sources.items()}

# 2. Calculate the sound level from each source at the new position
new_levels = {
    name: calculate_new_sound_level(data["Lp1"], new_distances[name])
    for name, data in sources.items()
}

# 3. Calculate the total combined sound level
total_level = calculate_total_sound_level(new_levels.values())

# --- Output ---
print(f"Your new position is {my_new_pos}.")
print("-" * 30)
for name, level in new_levels.items():
    print(f"The distance to the {name} is {new_distances[name]:.2f} meters.")
    print(f"The new sound level from the {name} is {level:.2f} dB.")
    print("-" * 30)

# Construct the final equation string for clarity
equation_parts = [f"10^({level:.2f}/10)" for level in new_levels.values()]
equation_string = " + ".join(equation_parts)
print("The final calculation is:")
print(f"Total dB = 10 * log10({equation_string})")
print(f"Total dB = {total_level:.2f}")

# Final Answer Wrapper
# <<<75.12>>>