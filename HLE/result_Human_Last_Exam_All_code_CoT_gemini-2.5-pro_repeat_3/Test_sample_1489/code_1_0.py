import re

def note_to_midi(note_name):
    """Converts a note name like 'A4' or 'C#5' to a MIDI number."""
    note_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    match = re.match(r'([A-G])(#?)(\d+)', note_name)
    note, sharp, octave = match.groups()
    
    midi_val = note_map[note]
    if sharp:
        midi_val += 1
    
    # MIDI formula: value + (octave + 1) * 12
    # C4 is middle C (60)
    midi_number = midi_val + (int(octave) + 1) * 12
    return midi_number

def midi_to_note(midi_number):
    """Converts a MIDI number to a note name."""
    note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    note_index = midi_number % 12
    octave = midi_number // 12 - 1
    return f"{note_names[note_index]}{octave}"

def main():
    # Initial sequence
    initial_notes = ["A4", "C5", "F4", "E4", "G4", "C4", "B4", "D4"]
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # 1. Transpose up Major Second (2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]

    # 2. Invert Around E4
    axis_e4 = note_to_midi("E4")
    midi_sequence = [(2 * axis_e4) - n for n in midi_sequence]

    # 3. Retrograde the sequence
    midi_sequence.reverse()

    # 4. Augment intervals by adding 3 semitones
    augmented_sequence = []
    if midi_sequence:
        augmented_sequence.append(midi_sequence[0])
        for i in range(len(midi_sequence) - 1):
            original_interval = midi_sequence[i+1] - midi_sequence[i]
            new_interval = original_interval + 3
            next_note = augmented_sequence[i] + new_interval
            augmented_sequence.append(next_note)
        midi_sequence = augmented_sequence

    # 5. Change to Dorian Mode Starting from D4
    # This is interpreted as transposing the melodic contour to a new starting note.
    # Calculate intervals from the first note of the current sequence.
    tonic = midi_sequence[0]
    intervals = [n - tonic for n in midi_sequence]
    # Apply these intervals to the new tonic, D4.
    new_tonic_d4 = note_to_midi("D4")
    midi_sequence = [new_tonic_d4 + i for i in intervals]
    
    # 6. Transpose down Minor Third (3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]

    # 7. Invert around F4
    axis_f4 = note_to_midi("F4")
    midi_sequence = [(2 * axis_f4) - n for n in midi_sequence]

    # 8. Transposed up one Octave (12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]

    # Convert final MIDI numbers back to note names
    final_notes = [midi_to_note(n) for n in midi_sequence]
    
    # Print the final sequence
    print(" ".join(final_notes))

main()