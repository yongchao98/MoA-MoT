def find_musical_joke_measures():
    """
    This function identifies and formats information about a famous musical joke
    in a piano concerto.

    The joke in question occurs in the third movement (Rondo) of Beethoven's
    Piano Concerto No. 5, Op. 73, the "Emperor" concerto. Towards the end,
    the cellos and basses enter with the main theme but are off by one beat.
    This "error" continues for two measures before they cleverly skip a beat
    to realign with the orchestra.

    - Composer: Ludwig van Beethoven
    - Opus Number: 73
    - Joke Start: The off-beat flourish begins in measure 317.
    - Joke Recovery: The cellos skip a beat in measure 320 to recover.
    """

    composer_surname = "Beethoven"
    opus_number = 73
    start_measure = 317
    end_measure = 320

    # Format the measure range as a hyphenated string
    measure_range = f"{start_measure}-{end_measure}"

    # Format the final answer as a comma-separated list
    # This fulfills the requirement to output each number (73, 317, 320)
    # in the final formatted response.
    final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

    print(final_answer)

find_musical_joke_measures()