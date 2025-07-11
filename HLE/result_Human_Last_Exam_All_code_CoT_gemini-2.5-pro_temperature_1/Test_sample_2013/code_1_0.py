def solve_musical_joke():
    """
    This function identifies and formats the details of a famous musical joke in a piano concerto.
    The joke in question is from Beethoven's Piano Concerto No. 5, "Emperor", at the transition
    to the third movement (Rondo). The cellos and basses anticipate the main theme by one beat,
    creating rhythmic dissonance.
    """
    composer_surname = "Beethoven"
    opus_number = 73
    start_measure = 98  # The measure the cellos/basses enter with the flourish
    end_measure = 100   # The measure the cellos/basses are silent to recover

    # The problem asks to output each number. We will format them into the final string.
    # Composer: Beethoven
    # Opus Number: 73
    # Measure Range: 98-100
    
    print(f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}")

solve_musical_joke()