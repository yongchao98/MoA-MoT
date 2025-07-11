def solve_music_riddle():
    """
    This function provides the answer to a music history riddle by identifying
    the composer, concerto, and specific measures of a famous musical joke.
    """
    
    # Step 1: Identify the composer.
    composer_surname = "Beethoven"
    
    # Step 2: Identify the opus number of the piano concerto.
    # The piece is the Piano Concerto No. 5 in E-flat major.
    opus_number = 73
    
    # Step 3: Identify the measure range of the joke in the third movement.
    # The cello/bass section enters early in measure 108.
    # They rest for measure 110 and recover with a correct entry in measure 111.
    start_measure = 108
    end_measure = 111
    
    # Step 4: Format the answer as a comma-separated string.
    # The format is: composer surname, opus number, range of measures.
    answer = f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}"
    
    print(answer)

solve_music_riddle()