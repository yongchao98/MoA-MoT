import re

def note_to_midi(note_str):
    """Converts a note string (e.g., 'A4', 'C#5') to a MIDI note number."""
    NOTE_MAP = {
        'C': 0, 'B#': 0,
        'C#': 1, 'Db': 1,
        'D': 2,
        'D#': 3, 'Eb': 3,
        'E': 4, 'Fb': 4,
        'F': 5, 'E#': 5,
        'F#': 6, 'Gb': 6,
        'G': 7,
        'G#': 8, 'Ab': 8,
        'A': 9,
        'A#': 10, 'Bb': 10,
        'B': 11, 'Cb': 11
    }
    # Use regex to robustly parse note name and octave
    match = re.match(r'([A-G][#b]?)(-?\d+)', note_str)
    if not match:
        raise ValueError(f"Invalid note format: {note_str}")
        
    name, octave_str = match.groups()
    octave = int(octave_str)
    
    # MIDI note number calculation where C4 is 60
    return NOTE_MAP[name] + (octave + 1) * 12

def midi_to_note(midi_num, use_sharps=True):
    """Converts a MIDI note number back to a note string."""
    NOTES_SHARP = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    NOTES_FLAT = ["C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"]
    
    notes = NOTES_SHARP if use_sharps else NOTES_FLAT
    
    octave = midi_num // 12 - 1
    note_index = midi_num % 12
    
    return notes[note_index] + str(octave)

def process_music_sequence():
    """
    Applies a series of eight transformations to a musical sequence.
    """
    # Initial sequence
    initial_sequence_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    midi_sequence = [note_to_midi(n) for n in initial_sequence_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]
    
    # 2. Invert Around E4
    axis_e4 = note_to_midi('E4')
    midi_sequence = [2 * axis_e4 - n for n in midi_sequence]

    # 3. Retrograde the sequence
    midi_sequence.reverse()
    
    # 4. Augment intervals by adding 3 semitones
    intervals = [midi_sequence[i] - midi_sequence[i-1] for i in range(1, len(midi_sequence))]
    augmented_intervals = [i + 3 for i in intervals]
    new_sequence = [midi_sequence[0]]
    for interval in augmented_intervals:
        new_sequence.append(new_sequence[-1] + interval)
    midi_sequence = new_sequence
    
    # 5. Change to Dorian Mode Starting from D4
    # (Interpreted as transposing the melodic contour to start on D4)
    intervals_from_root = [n - midi_sequence[0] for n in midi_sequence]
    new_root_d4 = note_to_midi('D4')
    midi_sequence = [new_root_d4 + i for i in intervals_from_root]

    # 6. Transpose down Minor Third (-3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]
    
    # 7. Invert around F4
    axis_f4 = note_to_midi('F4')
    midi_sequence = [2 * axis_f4 - n for n in midi_sequence]

    # 8. Transpose up one Octave (+12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]

    # Convert final MIDI sequence back to note names
    final_note_names = [midi_to_note(n) for n in midi_sequence]

    # Print the final result
    print(" ".join(final_note_names))

# Run the process and print the final sequence
process_music_sequence()