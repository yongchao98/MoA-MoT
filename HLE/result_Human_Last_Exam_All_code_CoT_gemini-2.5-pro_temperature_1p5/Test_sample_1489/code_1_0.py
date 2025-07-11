import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
    NOTES = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    # Handle two-character note names like 'C#'
    if len(note_name) > 1 and note_name[1] in ['#', 'b']:
        pass
    else: # Re-parse for single letter notes
        note_name = note_str[0]
        octave = int(note_str[1:])

    return NOTES[note_name] + (octave + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'A4'."""
    NOTES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return f"{NOTES[note_index]}{octave}"

def apply_transformations():
    """Applies the eight musical operations to the initial sequence."""
    
    # Initial sequence
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    notes = [note_to_midi(n) for n in initial_notes_str]

    # 1. Transpose up Major Second (2 semitones)
    notes = [n + 2 for n in notes]

    # 2. Invert Around E4
    pivot_e4 = note_to_midi('E4')
    notes = [(2 * pivot_e4) - n for n in notes]

    # 3. Retrograde the sequence
    notes = notes[::-1]

    # 4. Augment intervals by adding 3 semitones
    if len(notes) > 1:
        intervals = [notes[i] - notes[i-1] for i in range(1, len(notes))]
        augmented_intervals = [i + 3 for i in intervals]
        new_notes = [notes[0]]
        for interval in augmented_intervals:
            new_notes.append(new_notes[-1] + interval)
        notes = new_notes

    # 5. Change to Dorian Mode Starting from D4
    root_d4 = note_to_midi('D4')
    # Dorian scale relative intervals: W-H-W-W-W-H-W -> 0, 2, 3, 5, 7, 9, 10
    dorian_rel_intervals = {0, 2, 3, 5, 7, 9, 10}
    d_dorian_tones = {(root_d4 + i) % 12 for i in dorian_rel_intervals}

    new_notes = []
    for n in notes:
        rem = n % 12
        if rem in d_dorian_tones:
            new_notes.append(n)
        else:
            dist_down = min([(rem - t + 12) % 12 for t in d_dorian_tones])
            dist_up = min([(t - rem + 12) % 12 for t in d_dorian_tones])
            if dist_up < dist_down:
                new_notes.append(n + dist_up)
            else: # Rule: if equidistant, round down
                new_notes.append(n - dist_down)
    notes = new_notes

    # 6. Transpose down Minor Third (3 semitones)
    notes = [n - 3 for n in notes]

    # 7. Invert around F4
    pivot_f4 = note_to_midi('F4')
    notes = [(2 * pivot_f4) - n for n in notes]

    # 8. Transposed up one Octave (12 semitones)
    notes = [n + 12 for n in notes]

    # Convert final MIDI numbers back to note names
    final_note_names = [midi_to_note(n) for n in notes]
    
    print(" ".join(final_note_names))

# Run the transformations and print the result
apply_transformations()