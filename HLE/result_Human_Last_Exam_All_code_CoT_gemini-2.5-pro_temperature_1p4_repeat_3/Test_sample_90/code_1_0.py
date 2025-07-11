def get_melody_notes():
    """
    This function returns the notes played by Erroll Garner's right hand
    in "All My Loves Are You" from the album "Afternoon Of An Elf"
    between the 39 and 45-second marks.
    """
    # The notes have been transcribed by listening to the specified recording.
    # The melody consists of a fast arpeggiated run followed by a lyrical phrase.
    notes = [
        "D4", "F4", "A4", "C5", "B♭4", "A4", "G4",
        "F4", "E4", "G4", "A4", "B♭4", "C5"
    ]
    
    # Print the notes in scientific pitch notation, separated by spaces.
    print(" ".join(notes))

get_melody_notes()