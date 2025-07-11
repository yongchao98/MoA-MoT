import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a specific location.
    """
    # Step 1: Define initial positions and sound source data
    # My new position after walking 25m towards the construction site
    my_new_pos = (0, 25)

    # Sound sources: {name: {'pos': (x, y), 'spl_1m': dB at 1 meter}}
    sources = {
        'dog': {'pos': (-25, 0), 'spl_1m': 55},
        'train': {'pos': (50, 0), 'spl_1m': 110},
        'construction': {'pos': (0, 75), 'spl_1m': 90},
        'people': {'pos': (0, -10), 'spl_1m': 75}
    }

    new_spls = {}
    total_intensity = 0

    # Step 2 & 3: Calculate new SPL for each source and accumulate intensity
    for name, data in sources.items():
        source_pos = data['pos']
        spl_1m = data['spl_1m']

        # Calculate the distance from my new position to the sound source
        distance = math.sqrt((source_pos[0] - my_new_pos[0])**2 + (source_pos[1] - my_new_pos[1])**2)

        # Calculate the new sound pressure level (SPL) at this distance
        # Formula: SPL(r) = SPL(r_ref) - 20 * log10(r / r_ref) where r_ref=1
        spl_at_new_pos = spl_1m - 20 * math.log10(distance)
        new_spls[name] = spl_at_new_pos

        # Convert the new SPL to a value proportional to intensity
        # Formula: I = 10^(SPL / 10)
        intensity = 10**(spl_at_new_pos / 10)
        
        # Add to the total intensity
        total_intensity += intensity

    # Step 4: Convert total intensity back to a total sound level in dB
    # Formula: Total SPL = 10 * log10(Total Intensity)
    total_spl = 10 * math.log10(total_intensity)

    # Print the results and the final equation
    print("Sound level from each source at your new location:")
    print(f"Dog: {new_spls['dog']:.2f} dB")
    print(f"Train: {new_spls['train']:.2f} dB")
    print(f"Construction: {new_spls['construction']:.2f} dB")
    print(f"People: {new_spls['people']:.2f} dB")
    print("\n--------------------------------------------------")
    
    # Print the final equation with the calculated SPL values
    equation_str = (
        f"Total dB = 10 * log10("
        f"10^({new_spls['dog']:.2f}/10) + "
        f"10^({new_spls['train']:.2f}/10) + "
        f"10^({new_spls['construction']:.2f}/10) + "
        f"10^({new_spls['people']:.2f}/10))"
    )
    print("Final Equation:")
    print(equation_str)
    print("--------------------------------------------------")
    
    print(f"The total sound level you hear is: {total_spl:.2f} dB")
    
    # Final answer in the required format
    print(f"<<<{total_spl:.1f}>>>")

if __name__ == '__main__':
    calculate_total_sound_level()