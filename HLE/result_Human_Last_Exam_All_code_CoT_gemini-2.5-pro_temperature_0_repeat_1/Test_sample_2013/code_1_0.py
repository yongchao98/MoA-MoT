def solve_musical_puzzle():
    """
    This function identifies and formats the details of a famous musical joke
    in a piano concerto.
    """
    # The composer of the piano concerto.
    composer_surname = "Saint-Saëns"

    # The opus number of the piano concerto.
    opus_number = 22

    # The range of measures where the joke occurs.
    # The flourish starts in measure 103.
    # The cellos are off-beat for measures 104 and 105.
    # They skip a beat and recover in measure 106.
    measure_range = "103-106"

    # Assemble the final answer in the specified comma-separated format.
    final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

    print(final_answer)

solve_musical_puzzle()