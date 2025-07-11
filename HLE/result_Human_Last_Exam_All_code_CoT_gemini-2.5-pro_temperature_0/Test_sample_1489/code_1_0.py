import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
    note_names = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    if len(note_name) > 1 and note_name[1] == '#':
        pass
    elif len(note_name) > 1 and note_name[1] == 'b': # Handle flats just in case
        base_note = note_name[0]
        note_val = (note_names[base_note] - 1 + 12) % 12
    else:
        note_val = note_names[note_name]
        
    return (octave + 1) * 12 + note_val

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'A4'."""
    note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return note_names[note_index] + str(octave)

def snap_to_d_dorian(notes):
    """Snaps notes to the nearest note in the D Dorian scale."""
    # D Dorian MIDI remainders: D, E, F, G, A, B, C -> {2, 4, 5, 7, 9, 11, 0}
    dorian_remainders = {0, 2, 4, 5, 7, 9, 11}
    snapped_notes = []
    for note in notes:
        if note % 12 in dorian_remainders:
            snapped_notes.append(note)
            continue
        
        # Find the nearest Dorian note
        dist_down = 1
        while (note - dist_down) % 12 not in dorian_remainders:
            dist_down += 1
            
        dist_up = 1
        while (note + dist_up) % 12 not in dorian_remainders:
            dist_up += 1
        
        # If equidistant, round down (our convention)
        if dist_down <= dist_up:
            snapped_notes.append(note - dist_down)
        else:
            snapped_notes.append(note + dist_up)
            
    return snapped_notes

def main():
    # Initial sequence of notes
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    notes = [note_to_midi(n) for n in initial_notes_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    notes = [n + 2 for n in notes]

    # 2. Invert Around E4 (MIDI 64)
    pivot_e4 = note_to_midi('E4')
    notes = [(2 * pivot_e4) - n for n in notes]

    # 3. Retrograde the sequence
    notes.reverse()

    # 4. Augment the intervals by adding 3 semitones
    if len(notes) > 1:
        intervals = [notes[i] - notes[i-1] for i in range(1, len(notes))]
        augmented_intervals = [i + 3 for i in intervals]
        new_notes = [notes[0]]
        for interval in augmented_intervals:
            new_notes.append(new_notes[-1] + interval)
        notes = new_notes

    # 5. Change to Dorian Mode Starting from D4
    notes = snap_to_d_dorian(notes)

    # 6. Transpose down Minor Third (-3 semitones)
    notes = [n - 3 for n in notes]

    # 7. Invert around F4 (MIDI 65)
    pivot_f4 = note_to_midi('F4')
    notes = [(2 * pivot_f4) - n for n in notes]

    # 8. Transposed up one Octave (+12 semitones)
    notes = [n + 12 for n in notes]

    # Convert final MIDI numbers back to note names
    final_note_names = [midi_to_note(n) for n in notes]
    
    # Print the final result
    print(" ".join(final_note_names))

if __name__ == "__main__":
    main()