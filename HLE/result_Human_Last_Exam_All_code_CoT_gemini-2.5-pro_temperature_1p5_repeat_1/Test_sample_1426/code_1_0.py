def find_lowest_note():
    """
    This function identifies the lowest of the final three notes played by the 
    first violin in Dvořák's Symphony No. 2, Op. 4.

    The final bar for the first violins contains a single chord played 'divisi'
    (divided). The three notes that make up this chord are D5, F5, and B-flat5.
    To programmatically find the lowest note, we can assign them their MIDI
    pitch values, where a lower number corresponds to a lower pitch.
    """
    
    # The notes of the final chord and their MIDI values
    final_notes = {
        "D5": 74,
        "F5": 77,
        "B-flat5": 82
    }

    # Identify the note with the lowest MIDI value
    lowest_note = min(final_notes, key=final_notes.get)

    print("The final three notes for the first violin are in a single chord:")
    for note in final_notes:
        print(f"- {note}")
    
    print("\nOf these notes, the lowest one is:")
    print(lowest_note)

find_lowest_note()