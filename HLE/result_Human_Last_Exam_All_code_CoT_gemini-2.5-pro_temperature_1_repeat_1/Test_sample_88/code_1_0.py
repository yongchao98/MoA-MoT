def get_kaoru_improvisation():
    """
    This function returns the sequence of notes played by Kaoru Nishimi's
    right hand during his improvisation on "But Not For Me" in
    "Kids On the Slope" S1E4, measures 8-11.
    """
    # Measure 8: G5 E♭5 C5 B♭4
    # Measure 9: C6 B♭5 A♭5 G5
    # Measure 10: F5 E♭5 D5 C5
    # Measure 11: B♭4 C5 D5 E♭5
    notes = [
        "G5", "E♭5", "C5", "B♭4",
        "C6", "B♭5", "A♭5", "G5",
        "F5", "E♭5", "D5", "C5",
        "B♭4", "C5", "D5", "E♭5"
    ]
    print(" ".join(notes))

get_kaoru_improvisation()