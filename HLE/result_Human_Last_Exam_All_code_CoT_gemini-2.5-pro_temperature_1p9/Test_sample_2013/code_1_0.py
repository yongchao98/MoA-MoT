def find_musical_joke_measures():
    """
    Identifies the composer, opus, and measures for a famous musical joke,
    assuming the user meant the Haydn Cello Concerto instead of a piano concerto
    due to the specific details of the joke described.
    """

    # 1. Define the components of the answer
    composer_surname = "Haydn"
    opus_number = 101
    start_measure = 157
    end_measure = 159

    # 2. As requested, output each number that will appear in the final answer
    print(f"Opus Number: {opus_number}")
    print(f"Start Measure: {start_measure}")
    print(f"End Measure: {end_measure}")
    print("-" * 20)

    # 3. Format the final string as a comma-separated list
    measure_range = f"{start_measure}-{end_measure}"
    final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

    # 4. Print the final formatted answer
    print("The formatted answer is:")
    print(final_answer)

find_musical_joke_measures()
<<<Haydn, 101, 157-159>>>