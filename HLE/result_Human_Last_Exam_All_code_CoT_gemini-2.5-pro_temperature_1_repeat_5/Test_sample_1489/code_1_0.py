import math

def note_to_midi(note_str):
    """Converts a note name string (e.g., 'A4', 'C#5') to a MIDI number."""
    note_names = {"C": 0, "C#": 1, "D": 2, "D#": 3, "E": 4, "F": 5, "F#": 6, "G": 7, "G#": 8, "A": 9, "A#": 10, "B": 11}
    name = note_str[:-1]
    octave = int(note_str[-1])
    return note_names[name] + (octave + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number to a note name string."""
    note_names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return f"{note_names[note_index]}{octave}"

def apply_operations(initial_notes):
    """Applies the sequence of 8 musical operations."""

    # Convert initial note names to MIDI numbers
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # 1. Transpose up Major Second (+2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]

    # 2. Invert Around E4 (MIDI 64)
    axis_e4 = note_to_midi('E4')
    midi_sequence = [(2 * axis_e4) - n for n in midi_sequence]

    # 3. Retrograde the sequence
    midi_sequence.reverse()

    # 4. Augment intervals by adding 3 semitones
    if len(midi_sequence) > 1:
        intervals = [midi_sequence[i] - midi_sequence[i-1] for i in range(1, len(midi_sequence))]
        augmented_intervals = [i + 3 for i in intervals]
        new_sequence = [midi_sequence[0]]
        for interval in augmented_intervals:
            new_sequence.append(new_sequence[-1] + interval)
        midi_sequence = new_sequence

    # 5. Change to Dorian Mode Starting from D4
    d_dorian_pattern = [0, 2, 3, 5, 7, 9, 10] # Semitones from root
    d4_root = note_to_midi('D4')
    # Generate D Dorian scale notes over a wide range
    dorian_scale_notes = []
    for octave in range(-2, 10):
        for step in d_dorian_pattern:
            dorian_scale_notes.append(d4_root + octave * 12 + step - (d4_root % 12))

    # Sort the scale notes to be sure
    dorian_scale_notes.sort()
    
    new_sequence = []
    for note in midi_sequence:
        # Find the minimum distance to a scale note
        min_dist = min(abs(s - note) for s in dorian_scale_notes)
        # Find all scale notes at that minimum distance
        tied_notes = [s for s in dorian_scale_notes if abs(s - note) == min_dist]
        # Choose the highest note in case of a tie
        new_sequence.append(max(tied_notes))
    midi_sequence = new_sequence
    
    # 6. Transpose down Minor Third (-3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]
    
    # 7. Invert around F4 (MIDI 65)
    axis_f4 = note_to_midi('F4')
    midi_sequence = [(2 * axis_f4) - n for n in midi_sequence]

    # 8. Transposed up one Octave (+12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]
    
    # Convert final MIDI numbers back to note names
    final_notes = [midi_to_note(n) for n in midi_sequence]
    
    return final_notes

# Initial sequence of notes
start_sequence = ["A4", "C5", "F4", "E4", "G4", "C4", "B4", "D4"]

# Perform all operations
final_sequence = apply_operations(start_sequence)

# Print the final result as a space-separated string
print(" ".join(final_sequence))
