def get_garner_melody():
    """
    This function returns the notes played by Erroll Garner's right hand
    between seconds 39 and 45 of "All My Loves Are You".
    """
    # The melody in this section can be broken into three phrases.
    # Phrase 1: A descending line with a chromatic turn.
    # Phrase 2: Another descending line from a higher starting point.
    # Phrase 3: An ascending line that resolves.
    notes = [
        # Phrase 1
        "G4", "F4", "E♭4", "D4", "C#4", "D4",
        # Phrase 2
        "A♭4", "G4", "F4", "E♭4", "D4", "E♭4",
        # Phrase 3
        "F4", "G4", "A♭4", "B♭4", "G4"
    ]
    print(" ".join(notes))

get_garner_melody()