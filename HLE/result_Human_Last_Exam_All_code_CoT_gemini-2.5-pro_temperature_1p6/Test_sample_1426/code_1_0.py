def find_lowest_note():
    """
    This function finds the lowest of the last three notes played by the 
    first violin in Dvořák's Symphony No. 2.
    """
    
    # The last three notes are F5, D5, and B-flat4.
    # We represent them with their standard MIDI note numbers for comparison.
    # C4 (Middle C) is MIDI note 60.
    notes = {
        "F5": 77,
        "D5": 74,
        "B-flat4": 70
    }

    print("The final three notes for the first violin are F5, D5, and B-flat4.")
    print("To find the lowest one, we can compare their MIDI pitch values.")
    
    # "The final equation" showing each number
    print("\nThe equation to solve is finding the minimum of the following pitches:")
    for note, pitch in notes.items():
        print(f"Note: {note}, Pitch Value: {pitch}")

    # Find the key (note name) with the minimum value in the dictionary
    lowest_note_name = min(notes, key=notes.get)
    lowest_pitch_value = notes[lowest_note_name]

    print(f"\nComparing the pitch values {notes['F5']}, {notes['D5']}, and {notes['B-flat4']}, the lowest value is {lowest_pitch_value}.")
    print(f"\nTherefore, the lowest note is {lowest_note_name}.")

# Run the function
find_lowest_note()
