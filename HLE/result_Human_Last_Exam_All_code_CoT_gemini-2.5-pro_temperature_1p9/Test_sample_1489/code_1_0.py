import math

def note_to_midi(note_str):
    """Converts a note string (e.g., 'A4', 'C#5') to a MIDI number."""
    NOTES = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    if len(note_str) > 2 and note_str[1] in ['#', 'b']:
        name = note_str[:2]
        octave_str = note_str[2:]
    else:
        name = note_str[:1]
        octave_str = note_str[1:]
        
    octave = int(octave_str)
    
    return NOTES[name] + (octave + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number back to a note string."""
    NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    octave = midi_num // 12 - 1
    note_index = midi_num % 12
    return f"{NOTE_NAMES[note_index]}{octave}"

def is_white_key(midi_num):
    """Checks if a MIDI note corresponds to a white key (C Major/D Dorian scale)."""
    return (midi_num % 12) in [0, 2, 4, 5, 7, 9, 11]

def find_closest_dorian_note(midi_num):
    """Finds the closest note in the D Dorian scale (white keys). Ties are rounded up."""
    if is_white_key(midi_num):
        return midi_num
    
    lower = midi_num - 1
    while not is_white_key(lower):
        lower -= 1

    upper = midi_num + 1
    while not is_white_key(upper):
        upper += 1

    # Tie-breaking: round half up
    if (midi_num - lower) == (upper - midi_num):
        return upper
    elif (midi_num - lower) < (upper - midi_num):
        return lower
    else:
        return upper

def main():
    # Initial sequence of notes
    initial_notes = ['A4', 'C5', 'F4', 'E4', 'G4', 'C4', 'B4', 'D4']
    
    # Convert to MIDI numbers
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # Operation 1: Transpose up Major Second (+2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]
    
    # Operation 2: Invert Around E4 (MIDI 64)
    axis_e4 = note_to_midi('E4')
    midi_sequence = [(2 * axis_e4) - n for n in midi_sequence]

    # Operation 3: Retrograde the sequence
    midi_sequence = midi_sequence[::-1]
    
    # Operation 4: Augment the intervals between consecutive notes by adding 3 semitones
    augmented_sequence = []
    if midi_sequence:
        augmented_sequence.append(midi_sequence[0])
        for i in range(1, len(midi_sequence)):
            original_interval = midi_sequence[i] - midi_sequence[i-1]
            new_interval = original_interval + 3
            new_note = augmented_sequence[i-1] + new_interval
            augmented_sequence.append(new_note)
        midi_sequence = augmented_sequence

    # Operation 5: Change to Dorian Mode Starting from D4
    midi_sequence = [find_closest_dorian_note(n) for n in midi_sequence]

    # Operation 6: Transpose down Minor Third (-3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]

    # Operation 7: Invert around F4 (MIDI 65)
    axis_f4 = note_to_midi('F4')
    midi_sequence = [(2 * axis_f4) - n for n in midi_sequence]
    
    # Operation 8: Transposed up one Octave (+12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]
    
    # Convert final MIDI numbers back to note names
    final_notes = [midi_to_note(n) for n in midi_sequence]
    
    # Print the result
    print(" ".join(final_notes))

main()
<<<A5 D6 C#5 E5 B4 A4 C#5 F#4>>>