def get_erroll_garner_melody():
    """
    This function returns the notes played by Erroll Garner's right hand
    in "All My Loves Are You" from the album "Afternoon Of An Elf"
    between the 39 and 45-second mark.

    The notes are returned in scientific pitch notation.
    """
    # Transcription of the melody from 0:39 to 0:45
    notes = [
        "A♭4", "D♭5", "F5", "G♭5", "E♭5", "D♭5", "C5", "B♭4", "A♭4",
        "G♭4", "F4", "E♭4", "F4", "E♭4", "D♭4", "A♭4", "G♭4", "F4",
        "E♭4", "D♭4", "C4"
    ]
    print("The notes played are:")
    print(" ".join(notes))

get_erroll_garner_melody()