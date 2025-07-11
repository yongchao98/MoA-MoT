def get_kaoru_improvisation_notes():
    """
    Returns the sequence of notes played by Kaoru Nishimi's right hand
    in measures 8-11 of his "But Not For Me" improvisation.
    """
    # Notes transcribed from "Kids On the Slope" S1E4 (approx. 11:13)
    # Measure 8: E♭5 C5 A♭4 G4
    # Measure 9: F4 D4 E♭4 F4
    # Measure 10: G4 A♭4 A4 B♭4
    # Measure 11: C5 A♭4 G4 F4
    notes = [
        "E♭5", "C5", "A♭4", "G4",
        "F4", "D4", "E♭4", "F4",
        "G4", "A♭4", "A4", "B♭4",
        "C5", "A♭4", "G4", "F4"
    ]
    print(" ".join(notes))

get_kaoru_improvisation_notes()