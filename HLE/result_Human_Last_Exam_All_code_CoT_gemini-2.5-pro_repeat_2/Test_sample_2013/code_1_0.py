def find_musical_joke_measures():
    """
    Identifies and formats the details of a famous musical joke
    in a Beethoven piano concerto.
    """
    # Step 1: Define the composer's surname.
    composer_surname = "Beethoven"

    # Step 2: Define the opus number for Piano Concerto No. 5.
    opus_number = 73

    # Step 3: Define the start and end measures of the joke.
    # The cello/bass flourish begins in measure 350.
    # They are off-beat in measures 351 and 352.
    # They skip a beat and recover in measure 353.
    start_measure = 350
    end_measure = 353

    # Step 4: Format the measure range and the final answer string.
    # The f-string directly uses the number variables to build the final output.
    measure_range = f"{start_measure}-{end_measure}"
    final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

    # Step 5: Print the result.
    print(final_answer)

find_musical_joke_measures()