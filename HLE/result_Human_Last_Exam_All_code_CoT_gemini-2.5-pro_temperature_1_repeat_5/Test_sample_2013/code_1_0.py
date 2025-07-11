def find_musical_joke_measures():
    """
    This function identifies and prints the details of a famous musical joke
    in a piano concerto by Beethoven.
    """
    # Composer's surname
    composer = "Beethoven"
    
    # Opus number of the Piano Concerto No. 5 ("Emperor")
    opus_number = 73
    
    # The joke occurs in the third movement.
    # The cellos and basses begin their flourish one beat early in measure 342.
    # They remain off-beat for measures 343 and 344.
    # They skip a beat and recover in measure 345.
    start_measure = 342
    end_measure = 345
    
    # Format the range as a string
    measure_range = f"{start_measure}-{end_measure}"
    
    # Print the final formatted answer
    print(f"{composer}, {opus_number}, {measure_range}")

find_musical_joke_measures()