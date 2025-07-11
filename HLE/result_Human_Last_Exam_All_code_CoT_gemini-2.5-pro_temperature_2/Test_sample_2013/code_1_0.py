def find_musical_joke():
    """
    This function identifies and prints the details of a famous musical joke
    in a piano concerto by Beethoven.
    """
    # Details of the musical piece and the specific passage
    composer_surname = "Beethoven"
    opus_number = 15
    start_measure = 491
    end_measure = 493

    # Format the measure range as a string
    measure_range = f"{start_measure}-{end_measure}"

    # Print the final answer in the required format
    print(f"{composer_surname}, {opus_number}, {measure_range}")

find_musical_joke()
<<<Beethoven, 15, 491-493>>>