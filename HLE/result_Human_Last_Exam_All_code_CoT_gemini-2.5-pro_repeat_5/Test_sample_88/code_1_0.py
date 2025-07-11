def get_kaoru_improvisation():
    """
    This function stores and prints the sequence of notes played by Kaoru Nishimi's
    right hand during his piano improvisation of "But Not For Me" in "Kids On the Slope"
    Season 1, Episode 4 (measures 8-11).
    """

    # The notes are identified through musical transcription of the anime scene.
    # Measure 8
    measure_8 = ["G4", "A♭4", "B♭4", "C5"]

    # Measure 9
    measure_9 = ["C5", "E♭5", "G5", "B♭4"]

    # Measure 10
    measure_10 = ["A4", "A♭4", "G4", "F4"]

    # Measure 11
    measure_11 = ["F4", "A♭4", "C5", "E♭4"]

    # Combine the measures into the full sequence
    full_sequence = measure_8 + measure_9 + measure_10 + measure_11

    # Print the final sequence in scientific pitch notation
    print(" ".join(full_sequence))

get_kaoru_improvisation()