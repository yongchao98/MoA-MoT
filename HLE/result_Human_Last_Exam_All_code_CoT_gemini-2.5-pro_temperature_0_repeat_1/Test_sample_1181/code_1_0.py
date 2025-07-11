def find_solaris_character():
    """
    This function identifies the character from the 1972 movie 'Solaris'
    who is ashamed to miss the sound of rustling leaves on Earth.
    """
    # A dictionary mapping the multiple-choice options to the character names.
    characters = {
        'A': 'Kris',
        'B': 'Hari',
        'C': 'Snaut',
        'D': 'Sartorius',
        'E': 'Gibarian'
    }

    # In the film, Dr. Snaut expresses his weariness and nostalgia for Earth
    # to the main character, Kris Kelvin. He confesses his shame over
    # being homesick for simple, earthly things.
    correct_character_name = 'Snaut'

    # Find the correct option letter corresponding to the character's name.
    correct_option = None
    for option, name in characters.items():
        if name == correct_character_name:
            correct_option = option
            break

    if correct_option:
        print(f"The character who is ashamed to miss the sound of leaves rustling on Earth is Dr. {correct_character_name}.")
        print(f"This corresponds to answer choice: {correct_option}")
    else:
        print("The correct character was not found in the options.")

find_solaris_character()