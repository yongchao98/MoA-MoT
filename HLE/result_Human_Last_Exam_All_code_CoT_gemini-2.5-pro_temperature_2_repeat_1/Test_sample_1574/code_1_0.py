def find_character_info():
    """
    This function searches a small knowledge base to identify a character
    from Mikhail Bulgakov's "A Dog's Heart".
    """
    # A simple database of characters and their plot-relevant descriptions.
    character_database = {
        'Vasnetsova': 'A typist who Sharikov brings to the apartment and attempts to assault.',
        'Varvana': 'Not a character involved in this incident.',
        'Zina': 'The young maid of Professor Preobrazhensky, often harassed by Sharikov but not the subject of this specific assault incident.',
        'Maria': 'Not a prominent character in the novel.',
        'Darya': 'The cook for Professor Preobrazhensky.'
    }

    target_description = "attempts to assault"
    correct_character = None

    # Find the character matching the key plot point.
    for name, description in character_database.items():
        if target_description in description:
            correct_character = name
            break
    
    if correct_character:
        print(f"Based on the plot of 'A Dog's Heart', the woman Polygraf Sharikov attempted to assault was named Vasnetsova.")
        print(f"The correct option is: {correct_character}")
    else:
        print("Could not identify the character from the information.")

find_character_info()