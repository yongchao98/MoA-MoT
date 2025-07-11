def find_musical_joke_measures():
    """
    This function identifies and formats the details of a famous musical joke
    in a piano concerto where the cello section becomes desynchronized.
    """

    # The composer of the piano concerto.
    composer_surname = "Saint-Saëns"

    # The opus number of the concerto (Piano Concerto No. 2 in G minor).
    opus_number = 22

    # The measure where the cello's flourish begins, making it go long.
    start_measure = 67

    # The measure where the cello section skips a beat to get back in time.
    end_measure = 70

    # The joke spans from the start of the flourish to the recovery.
    measure_range = f"{start_measure}-{end_measure}"

    # Print the final answer in the required format.
    print(f"{composer_surname}, {opus_number}, {measure_range}")

find_musical_joke_measures()