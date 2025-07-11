def find_musical_joke_measures():
    """
    This function provides the details of a famous musical joke in a piano concerto.
    The information is based on established musicological knowledge.
    """
    
    # The composer of the piano concerto.
    composer_surname = "Saint-Saëns"
    
    # The opus number of the Piano Concerto No. 2 in G minor.
    opus_number = 22
    
    # The range of measures spanning the joke.
    # m. 97: Cellos play a 4-note flourish instead of 3, going a beat too long.
    # m. 98: Cellos are now off-beat.
    # m. 99: Cellos skip a beat to get back in sync with the orchestra.
    measure_range = "97-99"
    
    # Print the final formatted answer.
    print(f"{composer_surname}, {opus_number}, {measure_range}")

find_musical_joke_measures()