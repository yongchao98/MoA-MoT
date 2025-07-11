def get_kaoru_improvisation():
    """
    This function returns the sequence of notes from Kaoru Nishimi's
    improvisation of "But Not For Me" in "Kids On the Slope" S1E4.
    """
    # The notes for measures 8-11, inclusive.
    notes = [
        # Measure 8
        "G5", "A5", "G5",
        # Measure 9
        "F5", "E5", "D5", "C#5", "D5", "E5",
        # Measure 10
        "F5", "E5", "D5", "C#5", "D5", "C#5", "B4", "A4",
        # Measure 11
        "G#4", "A4", "B4", "C5"
    ]
    return " ".join(notes)

# Print the final sequence of notes.
print(get_kaoru_improvisation())