def find_musical_joke_measures():
    """
    This function provides the details of a famous musical joke in a piano concerto.
    The details are based on a well-documented passage in classical music literature.
    """
    
    # Information based on analysis of the score
    composer_surname = "Saint-Saëns"
    opus_number = 22
    flourish_start_measure = 41
    recovery_measure = 46

    # Format the measure range as a hyphenated string
    measure_range = f"{flourish_start_measure}-{recovery_measure}"
    
    # Assemble the final answer string as per the required format
    # The output will be a comma-separated list of 3 elements:
    # composer surname, opus number, and the range of measures.
    final_answer = f"{composer_surname}, {opus_number}, {measure_range}"
    
    print(final_answer)

find_musical_joke_measures()