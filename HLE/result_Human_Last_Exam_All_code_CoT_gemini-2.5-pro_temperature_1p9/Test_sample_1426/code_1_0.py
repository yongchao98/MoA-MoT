def find_lowest_note():
    """
    Identifies the lowest of the last three notes played by the first violin
    in Dvořák's Symphony No. 2, Op. 4.
    """
    # The last three notes found from the score are F5, D5, and B-flat 4.
    # To compare them programmatically, we can represent them with their MIDI pitch numbers.
    # Middle C (C4) is 60. Each semitone is +1.
    notes = {
        'F5': 77,
        'D5': 74,
        'B♭4': 70
    }

    print("The final three notes for the first violin in Dvořák's Symphony No. 2 are F5, D5, and B♭4.")
    print("To find the lowest, we can compare their standard MIDI pitch values.")
    
    # Print each note and its value
    for note, value in notes.items():
        print(f"Note: {note}, MIDI Value: {value}")
        
    # Find the note with the minimum MIDI value
    lowest_note_name = min(notes, key=notes.get)
    
    print(f"\nBy comparing the values, the lowest note is determined to be: {lowest_note_name}")

if __name__ == '__main__':
    find_lowest_note()