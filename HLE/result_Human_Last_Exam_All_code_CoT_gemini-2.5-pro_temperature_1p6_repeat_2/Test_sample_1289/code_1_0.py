def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """
    # Step 1 & 2: Define the note path in semitones from a root of C=0.
    # The scale is Maqam Zanjaran: Jins Hijaz on C + Jins Nahawand on F.
    # C=0, Db=1, E=4, F=5, G=7, Ab=8, Bb=10, C'=12.
    # The musician sings up the scale (C to C') and then descends from C' to F.
    # The octave C' is sung only once.
    sung_note_path_semitones = [
        0,  # C
        1,  # Db
        4,  # E
        5,  # F
        7,  # G
        8,  # Ab
        10, # Bb
        12, # C' (octave)
        10, # Bb
        8,  # Ab
        7,  # G
        5   # F
    ]

    # Step 3: Calculate the intervals.
    # An interval size is the difference in semitones / 2.
    # e.g., A semitone (1) becomes 0.5, a whole tone (2) becomes 1.
    intervals_in_tones = []
    for i in range(len(sung_note_path_semitones) - 1):
        note1 = sung_note_path_semitones[i]
        note2 = sung_note_path_semitones[i+1]
        semitone_difference = abs(note2 - note1)
        tone_value = semitone_difference / 2.0
        intervals_in_tones.append(tone_value)

    # Step 4: Format the output as a string: "{n1,n2,...}".
    # Numbers should be integers if they are whole, otherwise floats.
    formatted_intervals = []
    for val in intervals_in_tones:
        if val == int(val):
            formatted_intervals.append(str(int(val)))
        else:
            # The problem asks to round to the nearest quarter, but all our
            # intervals (0.5, 1.0, 1.5) are already exact multiples of 0.25.
            formatted_intervals.append(str(val))

    final_output = "{" + ",".join(formatted_intervals) + "}"
    print(final_output)

solve_maqam_intervals()
<<< {0.5,1.5,0.5,1,0.5,1,1,1,1,0.5,1} >>>