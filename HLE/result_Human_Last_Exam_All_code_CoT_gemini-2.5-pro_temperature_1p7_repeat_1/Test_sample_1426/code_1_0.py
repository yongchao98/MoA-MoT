def find_lowest_note():
    """
    This function identifies the lowest note among the last three notes played
    by the first violin in Dvořák's Symphony No. 2, Op. 4.
    """

    # Step 1: Define the last three notes based on score analysis.
    # The final passage for the first violin ends with the notes B5, B5, and B4.
    last_three_notes = ['B5', 'B5', 'B4']

    # Step 2: Prepare for programmatic comparison by converting notes to MIDI numbers.
    # This dictionary maps the 12 tones of the chromatic scale to a numeric value.
    NOTE_VALUES = {
        'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3, 'E': 4, 'F': 5,
        'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11
    }

    def note_to_midi(note_str):
        """Converts a note in scientific pitch notation (e.g., 'B4') to a MIDI number."""
        note_name = note_str[:-1]
        octave = int(note_str[-1])
        # MIDI formula: 12 * (octave + 1) + note_index
        return 12 * (octave + 1) + NOTE_VALUES[note_name]

    # Initialize variables to track the lowest note.
    lowest_note_name = None
    min_midi_value = float('inf')

    # Step 3: Print the notes being compared (the "equation").
    print("The final three notes played by the first violin are:")
    for note in last_three_notes:
        midi_value = note_to_midi(note)
        print(f"Note: {note} (MIDI value: {midi_value})")
        if midi_value < min_midi_value:
            min_midi_value = midi_value
            lowest_note_name = note

    # Step 4: Print the final answer.
    print(f"\nThe lowest of these three notes is {lowest_note_name}.")

find_lowest_note()