def get_erroll_garner_notes():
    """
    This function returns the notes played by Erroll Garner's right hand
    in "All My Loves Are You" between 39s and 45s.
    """
    # The notes identified from the recording in scientific pitch notation.
    # The passage consists of a descending line followed by a resolution and an ascending flourish.
    notes = [
        "D5", "C5", "A4", "F4",  # A descending line starting at ~39s
        "E5", "C5", "A4",         # A quick flourish at ~42s
        "D4",                     # A resolving note at ~43s
        "G4", "A4", "B4", "C5"     # An ascending run leading out of the phrase at ~44-45s
    ]

    # Print the notes separated by a space
    print(' '.join(notes))

get_erroll_garner_notes()