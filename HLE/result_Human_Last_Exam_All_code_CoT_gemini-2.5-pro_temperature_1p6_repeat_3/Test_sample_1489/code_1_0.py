import math

def note_to_midi(note_str):
    """Converts a note string like 'C#4' to a MIDI number."""
    note_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    base_note = note_name[0]
    midi_val = note_map[base_note]
    
    if len(note_name) > 1:
        if note_name[1] == '#':
            midi_val += 1
        elif note_name[1] == 'b':
            midi_val -= 1
            
    return midi_val + (octave + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'C#4'."""
    note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return f"{note_names[note_index]}{octave}"

def find_closest_in_scale(note, root, scale_degrees):
    """Finds the closest note in a given scale."""
    def is_in_scale(n):
        return (n - root) % 12 in scale_degrees

    if is_in_scale(note):
        return note

    offset = 1
    while True:
        # In case of a tie, this implementation will prefer the lower note.
        lower_note = note - offset
        if is_in_scale(lower_note):
            return lower_note
            
        upper_note = note + offset
        if is_in_scale(upper_note):
            return upper_note
            
        offset += 1

def main():
    # Initial sequence of eight musical notes
    initial_notes = ['A4', 'C5', 'F4', 'E4', 'G4', 'C4', 'B4', 'D4']
    
    # Convert to MIDI numbers for processing
    midi_seq = [note_to_midi(n) for n in initial_notes]

    # Operation 1: Transpose up Major Second (+2 semitones)
    midi_seq = [n + 2 for n in midi_seq]

    # Operation 2: Invert Around E4 (MIDI 64)
    inversion_axis_1 = note_to_midi('E4')
    midi_seq = [2 * inversion_axis_1 - n for n in midi_seq]

    # Operation 3: Retrograde the sequence
    midi_seq.reverse()

    # Operation 4: Augment the intervals between consecutive notes by adding 3 semitones
    augmented_seq = [midi_seq[0]]
    for i in range(1, len(midi_seq)):
        original_interval = midi_seq[i] - midi_seq[i-1]
        new_note = augmented_seq[-1] + original_interval + 3
        augmented_seq.append(new_note)
    midi_seq = augmented_seq

    # Operation 5: Change to Dorian Mode Starting from D4
    dorian_root = note_to_midi('D4')
    # T-S-T-T-T-S-T pattern in semitones from root: 0, 2, 3, 5, 7, 9, 10
    dorian_degrees = {0, 2, 3, 5, 7, 9, 10}
    midi_seq = [find_closest_in_scale(n, dorian_root, dorian_degrees) for n in midi_seq]

    # Operation 6: Transpose down Minor Third (-3 semitones)
    midi_seq = [n - 3 for n in midi_seq]

    # Operation 7: Invert around F4 (MIDI 65)
    inversion_axis_2 = note_to_midi('F4')
    midi_seq = [2 * inversion_axis_2 - n for n in midi_seq]
    
    # Operation 8: Transposed up one Octave (+12 semitones)
    midi_seq = [n + 12 for n in midi_seq]

    # Convert the final MIDI sequence back to note names
    final_notes = [midi_to_note(n) for n in midi_seq]

    # Print the final result as a space-separated string
    final_answer = " ".join(final_notes)
    print(final_answer)

if __name__ == "__main__":
    main()