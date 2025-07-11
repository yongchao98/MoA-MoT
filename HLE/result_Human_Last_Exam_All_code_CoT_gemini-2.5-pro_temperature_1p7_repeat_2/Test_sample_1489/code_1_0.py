import math

def note_to_midi(note_str):
    """Converts a note string (e.g., 'A4', 'C#5') to a MIDI number."""
    pitch_class_map = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    name = note_str.upper().replace('B#', 'C').replace('E#', 'F').replace('CB', 'B').replace('FB', 'E')
    
    if len(name) > 1 and (name[1] == '#' or name[1] == 'B'):
        pitch = name[:2]
        octave = int(name[2:])
    else:
        pitch = name[:1]
        octave = int(name[1:])
        
    return pitch_class_map[pitch] + 12 * (octave + 1)

def midi_to_note(midi_val):
    """Converts a MIDI number to a note string, using sharps for accidentals."""
    note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    pitch_class_index = midi_val % 12
    octave = midi_val // 12 - 1
    return f"{note_names[pitch_class_index]}{octave}"

def apply_operations(initial_notes_str):
    """Applies the full sequence of eight operations."""
    notes = [note_to_midi(n) for n in initial_notes_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    notes = [n + 2 for n in notes]

    # 2. Invert Around E4
    axis_e4 = note_to_midi('E4')
    notes = [2 * axis_e4 - n for n in notes]

    # 3. Retrograde the sequence
    notes.reverse()

    # 4. Augment intervals by adding 3 semitones
    notes = [notes[i] + 3 * i for i in range(len(notes))]

    # 5. Change to Dorian Mode Starting from D4
    dorian_intervals = [0, 2, 3, 5, 7, 9, 10]  # T-S-T-T-T-S-T
    root_d4 = note_to_midi('D4')
    dorian_pitch_classes = {(root_d4 + i) % 12 for i in dorian_intervals}
    
    scale_notes = sorted([i for i in range(128) if i % 12 in dorian_pitch_classes])
    
    new_notes = []
    import bisect
    for note in notes:
        idx = bisect.bisect_left(scale_notes, note)
        if idx == 0:
            nearest = scale_notes[0]
        elif idx == len(scale_notes):
            nearest = scale_notes[-1]
        else:
            before, after = scale_notes[idx - 1], scale_notes[idx]
            if note - before < after - note:
                nearest = before
            else: # In case of a tie, round up to the higher note
                nearest = after
        new_notes.append(nearest)
    notes = new_notes

    # 6. Transpose down Minor Third (-3 semitones)
    notes = [n - 3 for n in notes]

    # 7. Invert around F4
    axis_f4 = note_to_midi('F4')
    notes = [2 * axis_f4 - n for n in notes]
    
    # 8. Transpose up one Octave (+12 semitones)
    notes = [n + 12 for n in notes]

    # Convert final MIDI numbers to note names
    final_note_names = [midi_to_note(n) for n in notes]
    
    print(" ".join(final_note_names))

if __name__ == '__main__':
    initial_sequence = "A4 C5 F4 E4 G4 C4 B4 D4"
    apply_operations(initial_sequence)