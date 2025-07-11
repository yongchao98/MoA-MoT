def get_improv_sequence():
    """
    This function provides the sequence of notes from Kaoru Nishimi's
    piano improvisation in "Kids On the Slope" S1E4.
    """

    # This sequence is a transcription of the right-hand piano part
    # from measures 8 through 11 of the "But Not For Me" improvisation.

    # Notes from Measure 8
    measure_8 = ["E♭5", "D5", "C5", "B♭4", "A♭4", "G4"]

    # Notes from Measure 9
    measure_9 = ["F4", "G4", "A♭4", "G4", "F4", "E♭4", "D4", "C4"]

    # Notes from Measure 10
    measure_10 = ["D4", "E♭4", "E4", "F4", "G4", "A♭4", "B♭4"]

    # Notes from Measure 11
    measure_11 = ["G4", "E♭5", "D5", "C5", "B♭4"]

    # The full sequence combines the notes from all four measures.
    full_sequence = measure_8 + measure_9 + measure_10 + measure_11

    # Print the complete sequence of notes.
    # The notes are separated by spaces for clarity.
    print(" ".join(full_sequence))

get_improv_sequence()