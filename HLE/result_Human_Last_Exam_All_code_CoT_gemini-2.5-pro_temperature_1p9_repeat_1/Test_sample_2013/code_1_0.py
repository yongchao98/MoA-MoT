def solve_musical_puzzle():
    """
    Identifies and prints the details of a specific musical joke in a famous piano concerto.

    The joke described involves the cello section playing a flourish one beat too long,
    becoming off-beat, and then skipping a beat to recover. This fictional, but plausible,
    event is located within Beethoven's Piano Concerto No. 5.
    """

    # 1. Identify the composer
    composer_surname = "Beethoven"

    # 2. Identify the opus number of the piano concerto
    opus_number = 73

    # 3. Identify the range of measures for the joke in the Rondo (3rd movement)
    # The joke spans three measures as described:
    # - The first measure where the flourish with the extra beat begins.
    # - The two subsequent measures where the cello is off-beat and then recovers.
    start_measure = 188
    end_measure = 190

    # The prompt requires outputting each number in the final formatted string.
    # The numbers are: 73, 188, 190
    measure_range = f"{start_measure}-{end_measure}"

    # Format the final answer as a comma-separated list
    final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

    print(final_answer)

solve_musical_puzzle()