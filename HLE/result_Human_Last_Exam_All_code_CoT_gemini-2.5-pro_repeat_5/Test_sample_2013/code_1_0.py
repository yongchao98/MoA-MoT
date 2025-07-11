def solve_musical_joke():
    """
    This function identifies and prints the details of a famous musical joke
    in a piano concerto as per the user's request.
    """
    composer_surname = "Saint-Saëns"
    opus_number = 22
    start_measure = 368
    end_measure = 372

    # The final answer is a comma-separated list of three elements:
    # 1. Composer's surname
    # 2. Opus number
    # 3. The range of measures for the joke
    print(f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}")

solve_musical_joke()