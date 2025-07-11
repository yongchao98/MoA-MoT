import math

def main():
    """
    Calculates the total sound level from multiple sources at a new position.
    """
    # Step 1: Define initial positions and sound levels at 1 meter.
    # Your initial position is assumed to be (0, 0).
    # Directions: left is negative x, right is positive x, front is positive y, behind is negative y.
    sources = {
        'dog': {'pos': (-25, 0), 'Lp1': 55},
        'train': {'pos': (50, 0), 'Lp1': 110},
        'construction': {'pos': (0, 75), 'Lp1': 90},
        'people': {'pos': (0, -10), 'Lp1': 75}
    }

    # Your new position after walking 25 meters towards the construction at (0, 75)
    my_new_pos = (0, 25)

    # Step 2 & 3: Calculate the new distance and resulting sound level from each source.
    new_sound_levels = {}
    print("Calculating sound levels at the new location (25m towards construction):\n")
    for name, data in sources.items():
        # Calculate distance from the new position to the source
        distance = math.sqrt((data['pos'][0] - my_new_pos[0])**2 + (data['pos'][1] - my_new_pos[1])**2)
        
        # Calculate the sound pressure level (Lp) at that distance from the source
        # Lp2 = Lp1 - 20 * log10(d2/d1). Here d1 is 1m.
        new_lp = data['Lp1'] - 20 * math.log10(distance)
        new_sound_levels[name] = new_lp
        print(f"Distance to {name}: {distance:.2f} m, Sound Level: {new_lp:.2f} dB")

    # Step 4: Combine the sound levels
    # The formula is L_total = 10 * log10( Σ(10^(Li/10)) )
    sum_of_intensities = 0
    for lp in new_sound_levels.values():
        sum_of_intensities += 10**(lp / 10)
    
    total_lp = 10 * math.log10(sum_of_intensities)

    # Step 5: Print the final equation and the result
    lp_dog = new_sound_levels['dog']
    lp_train = new_sound_levels['train']
    lp_const = new_sound_levels['construction']
    lp_ppl = new_sound_levels['people']

    print("\n------------------------------------------------------------")
    print("The total sound level is found by combining the individual sound levels.")
    print("The final calculation is:")
    print(
        f"10 * log10(10^({lp_dog:.2f}/10) + 10^({lp_train:.2f}/10) + 10^({lp_const:.2f}/10) + 10^({lp_ppl:.2f}/10)) = {total_lp:.2f} dB"
    )
    print("------------------------------------------------------------")


if __name__ == "__main__":
    main()