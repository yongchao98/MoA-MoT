def solve_musical_riddle():
    """
    This function identifies and formats the answer to the musical riddle
    about the cello joke in a famous piano concerto.
    """
    # The composer of the Piano Concerto No. 5, "Emperor".
    composer_surname = "Beethoven"

    # The opus number for the Piano Concerto No. 5.
    opus_number = 73

    # The joke begins with the cellos' early flourish.
    start_measure = 351

    # The cellos skip a beat and recover, ending the "off-beat" passage.
    end_measure = 353

    # Format the final answer string as: composer, opus, start-end
    # The code explicitly prints each component of the final answer.
    print(f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}")

solve_musical_riddle()