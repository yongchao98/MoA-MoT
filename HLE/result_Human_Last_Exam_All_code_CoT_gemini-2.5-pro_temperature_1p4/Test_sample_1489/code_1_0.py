import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
    note_map = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    # Handle cases like 'C#4' or 'A#3'
    if len(note_name) > 1 and (note_name.endswith('#') or note_name.endswith('b')):
        pass # Note name is already correct
    else: # single letter note
        pass

    return 12 * (octave + 1) + note_map[note_name]

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'A4'."""
    note_map = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return f"{note_map[note_index]}{octave}"

def apply_operations(initial_notes):
    """Applies the sequence of 8 musical operations."""
    
    # --- Initial State ---
    # Sequence: A4 C5 F4 E4 G4 C4 B4 D4
    notes = [note_to_midi(n) for n in initial_notes]
    # MIDI: [69, 72, 65, 64, 67, 60, 71, 62]

    # --- 1. Transpose up Major Second (+2 semitones) ---
    notes = [n + 2 for n in notes]
    # Result: [71, 74, 67, 66, 69, 62, 73, 64]

    # --- 2. Invert Around E4 ---
    axis = note_to_midi('E4') # 64
    notes = [(2 * axis) - n for n in notes]
    # Result: [57, 54, 61, 62, 59, 66, 55, 64]

    # --- 3. Retrograde the sequence ---
    notes = notes[::-1]
    # Result: [64, 55, 66, 59, 62, 61, 54, 57]
    
    # --- 4. Augment intervals by +3 semitones ---
    intervals = [notes[i+1] - notes[i] for i in range(len(notes) - 1)]
    aug_intervals = [i + 3 for i in intervals]
    new_notes = [notes[0]]
    for interval in aug_intervals:
        new_notes.append(new_notes[-1] + interval)
    notes = new_notes
    # Result: [64, 58, 72, 68, 74, 76, 72, 78]

    # --- 5. Change to Dorian Mode Starting from D4 ---
    # D Dorian scale pitch classes (C, D, E, F, G, A, B): {0, 2, 4, 5, 7, 9, 11}
    dorian_pcs = {0, 2, 4, 5, 7, 9, 11}
    modal_notes = []
    for note in notes:
        pc = note % 12
        if pc in dorian_pcs:
            modal_notes.append(note)
        else:
            min_dist = 12
            best_note = note
            # Find the closest scale note, rounding down in case of a tie
            for scale_pc in dorian_pcs:
                dist_up = (scale_pc - pc + 12) % 12
                dist_down = (pc - scale_pc + 12) % 12
                if dist_down <= min_dist:
                    min_dist = dist_down
                    best_note = note - dist_down
                if dist_up < min_dist:
                    min_dist = dist_up
                    best_note = note + dist_up
            modal_notes.append(best_note)
    notes = modal_notes
    # Result: [64, 57, 72, 67, 74, 76, 72, 77]

    # --- 6. Transpose down Minor Third (-3 semitones) ---
    notes = [n - 3 for n in notes]
    # Result: [61, 54, 69, 64, 71, 73, 69, 74]

    # --- 7. Invert around F4 ---
    axis = note_to_midi('F4') # 65
    notes = [(2 * axis) - n for n in notes]
    # Result: [69, 76, 61, 66, 59, 57, 61, 56]

    # --- 8. Transposed up one Octave (+12 semitones) ---
    notes = [n + 12 for n in notes]
    # Result: [81, 88, 73, 78, 71, 69, 73, 68]
    
    # Convert final MIDI numbers back to note names
    final_note_names = [midi_to_note(n) for n in notes]
    return " ".join(final_note_names)

# Initial sequence of notes
initial_sequence = ["A4", "C5", "F4", "E4", "G4", "C4", "B4", "D4"]

# Calculate and print the final sequence
final_sequence = apply_operations(initial_sequence)
print(final_sequence)