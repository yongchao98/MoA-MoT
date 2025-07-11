def find_lowest_note():
    """
    Identifies the lowest of the last three notes in Dvořák's Symphony No. 2, 1st violin part.
    The score shows the final notes are a D major chord: D, F#, and A.
    """
    # A mapping of note names to their pitch class value (C=0, C#=1, etc.)
    note_to_pitch_class = {
        'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3,
        'E': 4, 'F': 5, 'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8,
        'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11
    }

    # The last three notes found in the score for the first violin.
    # The symphony ends on a D major chord.
    notes_in_final_chord = ["D5", "F#5", "A5"]

    print(f"The notes in the final chord are: {', '.join(notes_in_final_chord)}.")

    note_values = {}

    # Calculate a numerical pitch value for each note for comparison.
    # The formula for a MIDI-like value is: 12 * octave + pitch_class
    for note_str in notes_in_final_chord:
        # Separate the note letter(s) from the octave number
        if len(note_str) > 2 and note_str[1] in ['#', 'b']:
            name = note_str[:2]
            octave = int(note_str[2:])
        else:
            name = note_str[0]
            octave = int(note_str[1:])

        # The octave for MIDI calculation is octave + 1 in this notation system.
        pitch_value = 12 * (octave + 1) + note_to_pitch_class[name]
        note_values[note_str] = pitch_value
    
    print("To compare them, we can convert each note to a numerical pitch value:")
    for note, value in note_values.items():
        print(f"- {note} has a value of {value}")

    # Find the note with the lowest pitch value
    lowest_note = min(note_values, key=note_values.get)
    
    # Print the equation/comparison. The min() function compares each number.
    pitch_numbers = list(note_values.values())
    print(f"\nEquation to find the lowest note: min({pitch_numbers[0]}, {pitch_numbers[1]}, {pitch_numbers[2]})")
    print(f"The lowest pitch value is {note_values[lowest_note]}.")

    # The final answer is the note name itself, without the octave.
    lowest_note_name = lowest_note[0]
    if len(lowest_note) > 2 and lowest_note[1] in ['#', 'b']:
        lowest_note_name = lowest_note[:2]

    print(f"\nTherefore, the lowest note is: {lowest_note_name}")

find_lowest_note()