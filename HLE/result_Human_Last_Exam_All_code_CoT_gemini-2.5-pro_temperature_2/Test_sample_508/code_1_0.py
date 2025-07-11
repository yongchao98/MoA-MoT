import math

def calculate_total_sound_level():
    """
    Calculates the total sound level at a new location from multiple sources.
    """
    # 1. Define initial positions and sound levels (Lp1 at r1=1m)
    my_pos_initial = (0, 0)
    sources = {
        'dog': {'pos': (-25, 0), 'Lp1': 55},
        'train': {'pos': (50, 0), 'Lp1': 110},
        'construction': {'pos': (0, 75), 'Lp1': 90},
        'people': {'pos': (0, -10), 'Lp1': 75}
    }
    
    # 2. Determine the new position
    # Walk 25 meters towards the construction (positive y-direction)
    my_pos_new = (my_pos_initial[0], my_pos_initial[1] + 25)
    
    # Dictionary to hold the calculated sound levels at the new position
    new_sound_levels = {}
    
    print("Calculating sound levels at your new position (0, 25):")
    
    # 3. & 4. For each source, calculate the new distance and the new sound level (Lp2)
    for name, data in sources.items():
        source_pos = data['pos']
        Lp1 = data['Lp1']
        
        # Calculate new distance (r2) from new position to the source
        r2 = math.sqrt((source_pos[0] - my_pos_new[0])**2 + (source_pos[1] - my_pos_new[1])**2)
        
        # Calculate sound level (Lp2) at distance r2. r1 is 1 meter.
        # Lp2 = Lp1 - 20 * log10(r2 / r1) -> Lp2 = Lp1 - 20 * log10(r2)
        Lp2 = Lp1 - 20 * math.log10(r2)
        new_sound_levels[name] = Lp2
        
        print(f"- Distance to {name}: {r2:.2f} meters. Sound level: {Lp2:.2f} dB")
        
    # 5. Combine the sound levels
    # Lp_total = 10 * log10( SUM( 10^(Lpi/10) ) )
    
    Lp_dog = new_sound_levels['dog']
    Lp_train = new_sound_levels['train']
    Lp_construction = new_sound_levels['construction']
    Lp_people = new_sound_levels['people']

    # Sum of the intensities
    total_intensity_term = (
        10**(Lp_dog / 10) +
        10**(Lp_train / 10) +
        10**(Lp_construction / 10) +
        10**(Lp_people / 10)
    )
    
    total_db = 10 * math.log10(total_intensity_term)
    
    # 6. Present the final calculation and result
    print("\nTo find the total sound level, we use the formula:")
    print("Total_dB = 10 * log10(10^(dB1/10) + 10^(dB2/10) + ...)\n")
    print("Final Equation:")
    print(f"Total dB = 10 * log10(10^({Lp_dog:.2f}/10) + 10^({Lp_train:.2f}/10) + 10^({Lp_construction:.2f}/10) + 10^({Lp_people:.2f}/10))")
    
    print(f"\nThe total sound level you hear is: {total_db:.2f} dB")
    
    return total_db

# Run the calculation and get the final answer
final_answer = calculate_total_sound_level()
# The final answer will be printed within the function.
# This final line is for extracting the answer as requested.
# Do not remove it.
print(f"<<<{final_answer:.2f}>>>")