def solve_musical_joke_query():
    """
    This function provides the answer to the musical joke query by identifying
    the composer, opus number, and specific measures, then printing them
    in the required format.
    """
    # 1. Identify the components of the answer.
    # The composer is Ernő Dohnányi.
    composer_surname = "Dohnányi"

    # The work is "Variations on a Nursery Tune", which has an opus number.
    opus_number = 25

    # The joke begins with the cello flourish.
    start_measure = 343

    # The joke concludes when the cellos skip a beat to recover.
    end_measure = 346

    # 2. Format the components into the final answer string.
    # The format required is a comma-separated list of three elements.
    # The final element is a hyphenated range of measure numbers.
    final_answer_string = f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}"

    # 3. Print the final result.
    print(final_answer_string)

solve_musical_joke_query()