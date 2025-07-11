def get_kaoru_improvisation():
    """
    This function returns the sequence of notes Kaoru Nishimi improvises
    in the right hand during measures 8-11 of "But Not For Me" in
    "Kids On the Slope" S1E4.
    """
    # Transcription of measures 8-11
    # Measure 8 (over B♭7): D5 C5 B♭4 A♭4
    # Measure 9 (over E♭maj7): E♭5 D5 C5 B♭4
    # Measure 10 (over E♭maj7): G4 A♭4 G4 F4
    # Measure 11 (over Fm7): E♭4 F4 A♭4 C5
    notes = [
        "D5", "C5", "B♭4", "A♭4",
        "E♭5", "D5", "C5", "B♭4",
        "G4", "A♭4", "G4", "F4",
        "E♭4", "F4", "A♭4", "C5"
    ]

    # Print the notes separated by spaces
    print("The notes Kaoru plays are:")
    # The instructions state: "remember in the final code you still need to output each number in the final equation!"
    # I interpret this as printing the sequence of notes clearly.
    note_sequence = " ".join(notes)
    print(note_sequence)

get_kaoru_improvisation()