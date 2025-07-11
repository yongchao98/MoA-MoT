def find_musical_joke_info():
    """
    This function provides information about a specific musical joke in a piano concerto.
    """
    composer_surname = "Beethoven"
    
    # This refers to the Piano Concerto in E-flat major, WoO 4.
    # "WoO" stands for "Werk ohne Opuszahl" (Work without Opus number).
    # We will use the number '4' as its identifier.
    opus_number = 4
    
    # The range of measures spanning the joke in the Rondo movement.
    # m. 278: The cello begins its flourish.
    # mm. 279-280: The cello is rhythmically off-beat.
    # m. 281: The cello realigns with the orchestra.
    measure_range = "278-281"
    
    # Print the final formatted answer.
    # The final equation is the formatted string itself.
    print(f"{composer_surname}, {opus_number}, {measure_range}")

find_musical_joke_info()