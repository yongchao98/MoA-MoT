def find_lowest_note():
    """
    This script identifies the lowest note from the final chord played by the first violin
    in Dvořák's Symphony No. 2, Op. 4.
    """
    # Step 1: Define the musical notes from the score.
    # The final chord for the first violin consists of three notes: B-flat4, D5, and F5.
    notes = ["B-flat4", "D5", "F5"]

    # Step 2: Create a mapping of note names to a pitch value for comparison.
    # This helps in programmatically ordering the notes.
    pitch_map = {
        "C": 0, "C#": 1, "D-flat": 1, "D": 2, "D#": 3, "E-flat": 3,
        "E": 4, "F": 5, "F#": 6, "G-flat": 6, "G": 7, "G#": 8,
        "A-flat": 8, "A": 9, "A#": 10, "B-flat": 10, "B": 11
    }

    lowest_note = None
    min_value = float('inf')

    print("The final three notes for the first violin are:")
    
    # Step 3: Iterate through the notes, calculate their value, and find the minimum.
    for note in notes:
        # Output each note being processed.
        print(note)
        
        # Parse the note string to get pitch and octave.
        octave = int(note[-1])
        pitch_name = note[:-1]
        
        # Calculate a single comparable value for the note.
        value = octave * 12 + pitch_map[pitch_name]
        
        # Check if this note is the lowest found so far.
        if value < min_value:
            min_value = value
            lowest_note = note
            
    # Step 4: Print the lowest note.
    print("\nThe lowest of these three notes is:")
    print(lowest_note)

find_lowest_note()