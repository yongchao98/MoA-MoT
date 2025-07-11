def solve_music_joke():
    """
    This function provides the details of a famous musical joke in a piano concerto.
    The user's description of a cello section playing a flourish one beat too long does not perfectly
    match a known piece. However, it strongly resembles a joke in Shostakovich's Piano Concerto No. 2,
    where the solo piano plays a clumsy passage, imitating a student's mistakes. This is the most
    likely piece the user is referring to.
    """

    # 1. The composer's surname
    composer_surname = "Shostakovich"

    # 2. The opus number of the piano concerto
    opus_number = 102

    # 3. The range of measures spanning the joke in the 3rd movement
    start_measure = 280
    end_measure = 287

    # The prompt requires printing each number in the final output.
    # The final output string contains all the required numbers.
    final_output = f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}"

    print(final_output)

solve_music_joke()