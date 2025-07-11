import math

def calculate_total_sound_level():
    """
    Calculates the total sound level from multiple sources at a new listener position.
    """
    # 1. Initial setup
    # Listener's initial position
    my_initial_pos = {'x': 0, 'y': 0}

    # Sound sources: {name: {level_at_1m, position}}
    sources = {
        'dog':       {'L1': 55, 'pos': {'x': -25, 'y': 0}},
        'train':     {'L1': 110, 'pos': {'x': 50, 'y': 0}},
        'construction': {'L1': 90, 'pos': {'x': 0, 'y': 75}},
        'people':    {'L1': 75, 'pos': {'x': 0, 'y': -10}}
    }

    # 2. Determine listener's new position
    my_new_pos = {'x': my_initial_pos['x'], 'y': my_initial_pos['y'] + 25}

    print("Your new position is (x, y) = ({}, {}).".format(my_new_pos['x'], my_new_pos['y']))
    print("\nCalculating the sound level from each source at your new position:")

    # 3 & 4. Calculate new distances and new sound levels for each source
    new_levels_at_listener = {}
    sum_of_intensities = 0

    for name, data in sources.items():
        # Calculate new distance to the source
        distance = math.sqrt((data['pos']['x'] - my_new_pos['x'])**2 + (data['pos']['y'] - my_new_pos['y'])**2)
        
        # Calculate the sound level at the new distance
        # Formula: L2 = L1 - 20 * log10(d2/d1), where d1=1m
        # Using d1=1 simplifies the formula to L2 = L1 - 20 * log10(d2)
        L2 = data['L1'] - 20 * math.log10(distance)
        new_levels_at_listener[name] = L2
        
        # Add the relative intensity (10^(L2/10)) to the sum
        sum_of_intensities += 10**(L2 / 10)
        
        print(f"- The {name} is {distance:.2f} meters away, producing a sound level of {L2:.2f} dB.")

    # 5. Calculate total combined sound level
    total_db = 10 * math.log10(sum_of_intensities)

    # 6. Format and print the final output
    ldog = new_levels_at_listener['dog']
    ltrain = new_levels_at_listener['train']
    lconst = new_levels_at_listener['construction']
    lppl = new_levels_at_listener['people']

    print("\nTo find the total sound level, we sum the sound intensities and convert back to decibels.")
    final_equation = f"Total Sound Level = 10 * log10(10^({ldog:.2f}/10) + 10^({ltrain:.2f}/10) + 10^({lconst:.2f}/10) + 10^({lppl:.2f}/10)) = {total_db:.2f} dB"
    print("The final equation is:")
    print(final_equation)

calculate_total_sound_level()
<<<75.12>>>