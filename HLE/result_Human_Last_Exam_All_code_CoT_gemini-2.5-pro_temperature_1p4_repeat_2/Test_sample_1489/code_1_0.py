import math

def note_to_midi(note):
    """Converts a note string (e.g., 'A4', 'C#5') to a MIDI number."""
    note_map = {'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3, 'E': 4, 'F': 5, 'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11}
    name = note[:-1]
    octave = int(note[-1])
    return note_map[name] + (octave + 1) * 12

def midi_to_note(midi):
    """Converts a MIDI number to a note string using sharps."""
    note_map_sharp = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}
    octave = midi // 12 - 1
    note_index = midi % 12
    return f"{note_map_sharp[note_index]}{octave}"

def find_closest_in_scale(note_midi, scale_semitones):
    """Finds the closest MIDI note in a given scale. Tie-breaks by choosing the lower note."""
    min_dist = float('inf')
    closest_note = -1
    
    # Search for the closest note below and above the current note
    # Check lower
    for i in range(note_midi, note_midi - 12, -1):
        if (i % 12) in scale_semitones:
            dist = note_midi - i
            if dist < min_dist:
                min_dist = dist
                closest_note = i
            break
            
    # Check higher
    for i in range(note_midi + 1, note_midi + 12):
        if (i % 12) in scale_semitones:
            dist = i - note_midi
            # If distance is smaller, or equal (tie-break goes to lower note from previous check)
            if dist < min_dist:
                min_dist = dist
                closest_note = i
            break
            
    return closest_note

def main():
    # Initial sequence
    initial_notes = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # 1. Transpose up Major Second (2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]

    # 2. Invert Around E4
    axis_e4 = note_to_midi('E4')
    midi_sequence = [(2 * axis_e4) - n for n in midi_sequence]

    # 3. Retrograde the sequence
    midi_sequence.reverse()

    # 4. Augment intervals by adding 3 semitones
    if len(midi_sequence) > 1:
        augmented_sequence = [midi_sequence[0]]
        for i in range(len(midi_sequence) - 1):
            original_interval = midi_sequence[i+1] - midi_sequence[i]
            new_interval = original_interval + 3
            next_note = augmented_sequence[-1] + new_interval
            augmented_sequence.append(next_note)
        midi_sequence = augmented_sequence

    # 5. Change to Dorian Mode Starting from D4
    # D Dorian scale pattern: T-S-T-T-T-S-T. From D(2): D, E, F, G, A, B, C
    # Semitone values (relative to C): 2, 4, 5, 7, 9, 11, 0
    d_dorian_semitones = {0, 2, 4, 5, 7, 9, 11}
    midi_sequence = [find_closest_in_scale(n, d_dorian_semitones) for n in midi_sequence]

    # 6. Transpose down Minor Third (3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]

    # 7. Invert around F4
    axis_f4 = note_to_midi('F4')
    midi_sequence = [(2 * axis_f4) - n for n in midi_sequence]

    # 8. Transposed up one Octave (12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]

    # Convert final MIDI sequence back to notes and print
    final_notes = [midi_to_note(n) for n in midi_sequence]
    print(" ".join(final_notes))

if __name__ == "__main__":
    main()