import collections

def find_lowest_note():
    """
    Finds the lowest of the last three notes for the first violin in Dvorak's Symphony No. 2.
    """
    # Step 1: Define the last three notes and their corresponding MIDI values.
    # This data is obtained by reading the musical score.
    # The notes are F5, D4, and B-flat4.
    notes = {
        "D4": 62,
        "B-flat4": 70,
        "F5": 77
    }

    print("The last three notes for the first violin are D4, B-flat4, and F5.")
    print("To find the lowest, we can compare their numerical pitch values (MIDI numbers):")
    
    # Step 2: Print each note and its value.
    for note, value in notes.items():
        print(f"Note: {note}, Pitch Value: {value}")
        
    # Step 3: Find the note with the minimum MIDI value.
    lowest_note_name = min(notes, key=notes.get)
    
    # Step 4: Print the result.
    print(f"\nComparing the pitch values {list(notes.values())}, the lowest value is {notes[lowest_note_name]}.")
    print(f"Therefore, the lowest note is {lowest_note_name}.")

find_lowest_note()