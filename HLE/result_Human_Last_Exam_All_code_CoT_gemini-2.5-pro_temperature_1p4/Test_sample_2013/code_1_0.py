def find_musical_joke_details():
    """
    This function provides details about a famous musical joke in a piano concerto.
    The joke involves the cello section humorously going off-beat and then recovering.
    """
    # Details of the musical piece
    composer_surname = "Haydn"
    
    # The concerto is catalogued as Hob. XVIII:11. We will use the primary number.
    concerto_number = 11
    
    # The range of measures where the joke unfolds
    # Measure 131: The cello flourish begins and goes a beat too long.
    # Measures 132-133: The cellos are off-beat.
    # Measure 134: The cellos skip a beat to recover.
    start_measure = 131
    end_measure = 134
    
    measure_range = f"{start_measure}-{end_measure}"
    
    # Print the formatted answer
    print(f"{composer_surname}, {concerto_number}, {measure_range}")

find_musical_joke_details()