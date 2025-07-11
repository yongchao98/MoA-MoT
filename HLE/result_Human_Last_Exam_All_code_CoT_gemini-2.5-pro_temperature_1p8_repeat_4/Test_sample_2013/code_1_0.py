def solve_music_riddle():
    """
    This function identifies a famous musical joke in a piano concerto
    and prints the composer, opus, and measure range.
    """
    composer = "Beethoven"
    opus_number = 73
    start_measure = 99
    end_measure = 102

    # Format the final answer as a comma-separated string
    # Composer Surname, Opus Number, Measure Range
    result = f"{composer}, {opus_number}, {start_measure}-{end_measure}"
    print(result)

solve_music_riddle()