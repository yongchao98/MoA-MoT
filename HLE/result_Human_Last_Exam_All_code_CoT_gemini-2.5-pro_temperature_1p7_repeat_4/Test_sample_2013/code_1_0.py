def find_musical_joke_measures():
    """
    This function provides information about a famous musical joke
    in a piano concerto by Beethoven.
    """
    
    # Composer's surname
    composer = "Beethoven"
    
    # Opus number for the Piano Concerto No. 1 in C Major
    opus_number = 15
    
    # The joke begins in the 3rd movement (Rondo).
    # The cello flourish starts in measure 401.
    start_measure = 401
    
    # The cellos are off-beat through measures 402 and 403.
    # They skip a beat in measure 404 to recover.
    # The range is therefore from the start of the flourish to the recovery.
    end_measure = 404
    
    # Format the final string as per the user's request.
    # The format is: composer surname, opus number, start_measure-end_measure
    result = f"{composer}, {opus_number}, {start_measure}-{end_measure}"
    
    print(result)

find_musical_joke_measures()
<<<Beethoven, 15, 401-404>>>