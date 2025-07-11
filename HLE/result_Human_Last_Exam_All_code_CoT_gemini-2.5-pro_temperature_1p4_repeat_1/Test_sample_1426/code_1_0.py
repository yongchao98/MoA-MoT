import re

def get_pitch_value(note_name):
    """Converts a musical note name (e.g., 'F#5') into a numerical pitch value."""
    # Mapping of note names to their position in an octave (0-11)
    note_map = {
        'C': 0, 'C#': 1, 'DB': 1,
        'D': 2, 'D#': 3, 'EB': 3,
        'E': 4, 'FB': 4,
        'F': 5, 'E#': 5, 'F#': 6, 'GB': 6,
        'G': 7, 'G#': 8, 'AB': 8,
        'A': 9, 'A#': 10, 'BB': 10,
        'B': 11, 'CB': 11
    }
    
    # Use regex to separate the note letter(s) and the octave number
    match = re.match(r'([A-G][#b]?)(\d+)', note_name)
    if not match:
        raise ValueError(f"Invalid note format: {note_name}")
        
    name, octave = match.groups()
    octave = int(octave)
    
    # Calculate the pitch value. This formula is similar to MIDI note calculation.
    return note_map[name.upper()] + octave * 12

def find_lowest_note():
    """
    Identifies the lowest note from the final chord in Dvořák's Symphony No. 2 for the first violin.
    """
    # The first violin part ends on a three-note chord.
    final_notes = ["B3", "F#5", "B5"]

    # Find the note with the minimum pitch value
    lowest_note = min(final_notes, key=get_pitch_value)

    print(f"The first violin ends Dvořák's Symphony No. 2 on a chord with three notes.")
    print(f"The notes in the final chord are: {', '.join(final_notes)}")
    print(f"Of these, the lowest note is {lowest_note}.")

if __name__ == "__main__":
    find_lowest_note()
