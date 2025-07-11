def find_musical_joke():
    """
    This function identifies and formats information about a famous musical joke
    in a piano concerto by Beethoven.
    """
    # Step 1: Define the components of the answer.
    # The piece is Beethoven's Piano Concerto No. 5.
    composer_surname = "Beethoven"
    
    # The opus number for this concerto is 73.
    opus_number = 73
    
    # Step 2: Define the measure range of the musical joke.
    # In the third movement (Rondo), the cello and bass flourish begins at measure 284.
    # They enter a beat late at measure 288 and remain off-beat through measure 289.
    # They recover in measure 290 by shortening a note.
    # The range is therefore from the start of the flourish to the measure of recovery.
    start_measure = 284
    end_measure = 290
    
    # Step 3: Format and print the final answer string.
    # The format required is: composer, opus, start-end
    print(f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}")

find_musical_joke()