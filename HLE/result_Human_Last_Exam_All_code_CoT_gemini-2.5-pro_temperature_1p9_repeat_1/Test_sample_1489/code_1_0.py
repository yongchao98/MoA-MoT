import math

def solve_music_sequence():
    """
    Applies a series of eight musical operations to a sequence of notes
    and prints the final result.
    """

    # --- Helper Functions ---

    def note_to_midi(note_str):
        """Converts a note string like 'A#4' or 'C5' to a MIDI number."""
        note_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
        name = note_str[0].upper()
        octave_part = note_str[1:]
        
        accidental_val = 0
        if len(octave_part) > 1 and octave_part[0] in '#b':
            if octave_part[0] == '#':
                accidental_val = 1
            else:  # 'b'
                accidental_val = -1
            octave = int(octave_part[1:])
        else:
            octave = int(octave_part)
            
        base_val = note_map[name]
        return (octave + 1) * 12 + base_val + accidental_val

    def midi_to_note(midi_num):
        """Converts a MIDI number back to a note string."""
        note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
        octave = (midi_num // 12) - 1
        note_index = midi_num % 12
        note_name = note_names[note_index]
        return f"{note_name}{octave}"

    # --- Initial State ---
    
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    notes = initial_notes_str.split()
    midi_sequence = [note_to_midi(n) for n in notes]
    # Initial sequence: [69, 72, 65, 64, 67, 60, 71, 62]

    # --- Applying 8 Operations ---

    # 1. Transpose up Major Second (+2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]
    # Result: [71, 74, 67, 66, 69, 62, 73, 64]

    # 2. Invert Around E4 (MIDI 64)
    axis_e4 = note_to_midi('E4')
    midi_sequence = [2 * axis_e4 - n for n in midi_sequence]
    # Result: [57, 54, 61, 62, 59, 66, 55, 64]

    # 3. Retrograde the sequence
    midi_sequence.reverse()
    # Result: [64, 55, 66, 59, 62, 61, 54, 57]
    
    # 4. Augment intervals by adding 3 semitones
    intervals = [midi_sequence[i+1] - midi_sequence[i] for i in range(len(midi_sequence) - 1)]
    augmented_intervals = [i + 3 for i in intervals]
    new_sequence = [midi_sequence[0]]
    for interval in augmented_intervals:
        new_sequence.append(new_sequence[-1] + interval)
    midi_sequence = new_sequence
    # Result: [64, 58, 72, 68, 74, 76, 72, 78]

    # 5. Change to Dorian Mode Starting from D4
    d_dorian_scale_notes = set()
    d4_midi = note_to_midi('D4')
    pattern = [2, 1, 2, 2, 2, 1, 2] # Dorian scale steps (W-H-W-W-W-H-W)
    # Generate scale notes in a wide range around the sequence's pitches
    current_note = d4_midi - 24 
    for i in range(60):
        d_dorian_scale_notes.add(current_note)
        current_note += pattern[i % 7]

    quantized_sequence = []
    for note in midi_sequence:
        if note in d_dorian_scale_notes:
            quantized_sequence.append(note)
        else:
            # Find nearest note(s) in scale
            min_diff = min(abs(note - scale_note) for scale_note in d_dorian_scale_notes)
            # Get all candidates with the minimum difference
            candidates = [sn for sn in d_dorian_scale_notes if abs(note - sn) == min_diff]
            # Tie-breaking rule: choose the lower note
            quantized_sequence.append(min(candidates))
    midi_sequence = quantized_sequence
    # Result: [64, 57, 72, 67, 74, 76, 72, 77]
    
    # 6. Transpose down Minor Third (-3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]
    # Result: [61, 54, 69, 64, 71, 73, 69, 74]

    # 7. Invert around F4 (MIDI 65)
    axis_f4 = note_to_midi('F4')
    midi_sequence = [2 * axis_f4 - n for n in midi_sequence]
    # Result: [69, 76, 61, 66, 59, 57, 61, 56]
    
    # 8. Transpose up one Octave (+12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]
    # Result: [81, 88, 73, 78, 71, 69, 73, 68]

    # --- Final Conversion and Output ---
    final_notes = [midi_to_note(n) for n in midi_sequence]
    print(" ".join(final_notes))

solve_music_sequence()