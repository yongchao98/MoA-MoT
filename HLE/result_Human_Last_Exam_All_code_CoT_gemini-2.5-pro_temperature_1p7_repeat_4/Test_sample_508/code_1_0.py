import math

def calculate_total_sound_level():
    """
    This function calculates the total sound level from multiple sources at a specific location.
    """
    # 1. Define initial positions and sound levels of sources (L1 at d1=1m).
    # We set the initial observer position at the origin (0, 0).
    sources = {
        'dog':          {'L1': 55, 'pos': (-25, 0)},
        'train':        {'L1': 110, 'pos': (50, 0)},
        'construction': {'L1': 90, 'pos': (0, 75)},
        'people':       {'L1': 75, 'pos': (0, -10)}
    }
    d1 = 1.0  # Initial distance for L1 measurements is 1 meter.

    # Your new position after walking 25m towards the construction site.
    my_new_pos = (0, 25)

    sound_levels_at_new_pos = {}
    relative_intensities = {}

    print("Step 1: Calculate the sound level (dB) from each source at your new location.")
    # 2. For each source, calculate the new distance (d2) and the new sound level (L2).
    for name, data in sources.items():
        source_pos = data['pos']
        L1 = data['L1']
        
        # Calculate new distance d2
        d2 = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)
        
        # Calculate new sound level L2
        L2 = L1 - 20 * math.log10(d2 / d1)
        
        sound_levels_at_new_pos[name] = L2
        print(f"- The {name}'s sound at {d2:.2f} meters is {L2:.2f} dB.")

    # 3. Convert each L2 to a relative intensity, sum them, and convert back to a total dB level.
    for name, L2 in sound_levels_at_new_pos.items():
        intensity = 10**(L2 / 10)
        relative_intensities[name] = intensity

    total_intensity = sum(relative_intensities.values())
    total_db = 10 * math.log10(total_intensity)

    # Unpack values for the final equation printout
    L2_dog = sound_levels_at_new_pos['dog']
    L2_train = sound_levels_at_new_pos['train']
    L2_construction = sound_levels_at_new_pos['construction']
    L2_people = sound_levels_at_new_pos['people']
    
    print("\nStep 2: Combine the sound levels.")
    print("This is done by converting dB to intensity, summing intensities, and converting back to dB.")
    print("\nThe final calculation is:")

    # Using f-strings to format the output with the calculated values
    equation_string = (
        f"Total dB = 10 * log10("
        f"10^({L2_dog:.2f}/10) + "
        f"10^({L2_train:.2f}/10) + "
        f"10^({L2_construction:.2f}/10) + "
        f"10^({L2_people:.2f}/10))"
    )
    print(f"{equation_string} = {total_db:.2f} dB")


calculate_total_sound_level()
<<<75.14>>>