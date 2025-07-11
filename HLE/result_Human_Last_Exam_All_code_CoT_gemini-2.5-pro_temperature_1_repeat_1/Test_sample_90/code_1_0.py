def get_erroll_garner_notes():
    """
    This function returns the notes played by Erroll Garner's right hand
    in "All My Loves Are You" between 0:39 and 0:45.
    """
    # The notes were determined by listening to the specified recording.
    # The sequence is a descending line followed by a brief ascending turnaround and run.
    notes = [
        "G5", "F5", "E5", "D5",
        "C5", "B♭4", "A4", "G4",
        "A4", "B♭4", "C5",
        "A4", "B♭4", "C5", "D5"
    ]
    print(" ".join(notes))

get_erroll_garner_notes()