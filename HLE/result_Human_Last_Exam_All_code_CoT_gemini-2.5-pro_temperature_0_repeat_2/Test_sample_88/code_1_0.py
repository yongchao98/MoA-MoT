def print_piano_improvisation():
    """
    Prints the sequence of notes from Kaoru Nishimi's piano improvisation
    on "But Not For Me" (measures 8-11) from "Kids On the Slope".
    """
    # The sequence of notes for measures 8-11, inclusive.
    # Measure 8: G4 B♭4 C5 E♭5
    # Measure 9: F5 E♭5 D5 C5 B♭4 A♭4 G4 F4
    # Measure 10: A♭4 B♭4 C5 D♭5 D5 E♭5 F5 G5
    # Measure 11: F5 E♭5 D5 C5
    notes = [
        "G4", "B♭4", "C5", "E♭5",
        "F5", "E♭5", "D5", "C5", "B♭4", "A♭4", "G4", "F4",
        "A♭4", "B♭4", "C5", "D♭5", "D5", "E♭5", "F5", "G5",
        "F5", "E♭5", "D5", "C5"
    ]

    # Print each note in the sequence, separated by a space.
    print(" ".join(notes))

print_piano_improvisation()