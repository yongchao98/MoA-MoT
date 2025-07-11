def print_erroll_garner_melody():
    """
    Prints the transcribed notes from Erroll Garner's "All My Loves Are You"
    (from the album "Afternoon Of An Elf") between seconds 39 and 45.
    """
    notes = [
        # Ascending run starting at ~39s
        "D5", "E5", "F5", "G5",
        # Held note and descending fill at ~40s
        "A5", "G5", "F5",
        # Ascending run at ~41s
        "C5", "D5", "E♭5", "F5",
        # Held note and descending phrase at ~42s
        "G5", "F5", "E♭5", "D5",
        # Lower phrase at ~43s
        "C5", "B♭4", "A4", "G4", "A4", "B♭4",
        # Resolution note at ~44s
        "C5"
    ]
    
    # Print the notes in a single line, separated by spaces
    print(" ".join(notes))

print_erroll_garner_melody()