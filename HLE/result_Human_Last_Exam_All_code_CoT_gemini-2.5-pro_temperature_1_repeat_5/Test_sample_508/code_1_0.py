import math

def calculate_total_sound_level():
    """
    This script calculates the total sound level from multiple sources at a specific location.
    """

    # Step 1 & 2: Define sound sources and the observer's new position.
    my_new_pos = (0, 25)

    sources = {
        'dog':          {'Lp1': 55,  'pos': (-25, 0)},
        'train':        {'Lp1': 110, 'pos': (50, 0)},
        'construction': {'Lp1': 90,  'pos': (0, 75)},
        'people':       {'Lp1': 75,  'pos': (0, -10)}
    }

    new_Lp_values = {}
    total_intensity = 0

    # Steps 3, 4, 5: Calculate distance, new Lp, and intensity for each source.
    for name, data in sources.items():
        # Calculate new distance from the observer's new position.
        distance = math.sqrt((data['pos'][0] - my_new_pos[0])**2 + (data['pos'][1] - my_new_pos[1])**2)

        # Calculate the attenuated sound level (Lp) at the new distance.
        # Formula: Lp2 = Lp1 - 20 * log10(distance) since reference distance is 1m.
        Lp_at_new_distance = data['Lp1'] - 20 * math.log10(distance)
        new_Lp_values[name] = Lp_at_new_distance

        # Convert the new sound level to relative intensity and add to the total.
        intensity = 10**(Lp_at_new_distance / 10)
        total_intensity += intensity

    # Step 6: Convert the total intensity back to a total sound level in dB.
    total_Lp = 10 * math.log10(total_intensity)

    # Print the final equation with all the calculated values and the final result.
    lp_dog = new_Lp_values['dog']
    lp_train = new_Lp_values['train']
    lp_construction = new_Lp_values['construction']
    lp_people = new_Lp_values['people']
    
    # Final output showing the equation with calculated numbers.
    print(
        f"The total sound level is calculated by converting each source's dB level at your new location to intensity, summing them, and converting back to dB.\n"
        f"\nFinal Equation:\n"
        f"Total dB = 10 * log10(10^({lp_dog:.2f}/10) + 10^({lp_train:.2f}/10) + 10^({lp_construction:.2f}/10) + 10^({lp_people:.2f}/10)) = {total_Lp:.2f} dB"
    )
    
    # Return the final value for the "answer" block.
    return total_Lp

# Execute the function and capture the final result
final_answer = calculate_total_sound_level()
# The final result is printed inside the function, this is just to capture it for the special format.
# <<<final_answer>>>