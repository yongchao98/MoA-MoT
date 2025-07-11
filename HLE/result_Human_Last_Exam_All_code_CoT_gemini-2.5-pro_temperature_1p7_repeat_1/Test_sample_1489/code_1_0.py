import math

def note_to_midi(note_str):
    """Converts a note string (e.g., 'A4', 'C#5') to a MIDI number."""
    note_map = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    name = note_str[:-1].upper()
    octave = int(note_str[-1])
    return note_map[name] + (octave + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number back to a note string."""
    note_map = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return f"{note_map[note_index]}{octave}"

def main():
    # Initial sequence of notes
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    
    # --- Step 0: Convert to MIDI ---
    notes = [note_to_midi(n) for n in initial_notes_str]
    # Sequence: [69, 72, 65, 64, 67, 60, 71, 62]

    # --- Step 1: Transpose up Major Second (+2 semitones) ---
    notes = [n + 2 for n in notes]
    # Result: [71, 74, 67, 66, 69, 62, 73, 64]

    # --- Step 2: Invert Around E4 (MIDI 64) ---
    axis_e4 = note_to_midi('E4') # 64
    notes = [2 * axis_e4 - n for n in notes]
    # Result: [57, 54, 61, 62, 59, 66, 55, 64]

    # --- Step 3: Retrograde the sequence ---
    notes.reverse()
    # Result: [64, 55, 66, 59, 62, 61, 54, 57]

    # --- Step 4: Augment intervals by +3 semitones ---
    intervals = [notes[i+1] - notes[i] for i in range(len(notes) - 1)]
    augmented_intervals = [i + 3 for i in intervals]
    new_notes = [notes[0]]
    for interval in augmented_intervals:
        new_notes.append(new_notes[-1] + interval)
    notes = new_notes
    # Result: [64, 58, 72, 68, 74, 76, 72, 78]

    # --- Step 5: Change to Dorian Mode Starting from D4 ---
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10} # Semitones from root
    d4_root = note_to_midi('D4') # 62
    
    # Generate a wide range of D Dorian notes
    dorian_scale_notes = set()
    for i in range(128):
        if (i - d4_root) % 12 in dorian_intervals:
            dorian_scale_notes.add(i)
    
    dorian_notes_sorted = sorted(list(dorian_scale_notes))
    
    mapped_notes = []
    for note in notes:
        if note in dorian_scale_notes:
            mapped_notes.append(note)
        else:
            # Find the closest note in the Dorian scale
            min_dist = float('inf')
            closest_note = -1
            for scale_note in dorian_notes_sorted:
                dist = abs(note - scale_note)
                if dist < min_dist:
                    min_dist = dist
                    closest_note = scale_note
                # Tie-breaking rule: choose the lower note
                elif dist == min_dist:
                    closest_note = min(closest_note, scale_note)
            mapped_notes.append(closest_note)
    notes = mapped_notes
    # Result: [64, 57, 72, 67, 74, 76, 72, 77]

    # --- Step 6: Transpose down Minor Third (-3 semitones) ---
    notes = [n - 3 for n in notes]
    # Result: [61, 54, 69, 64, 71, 73, 69, 74]

    # --- Step 7: Invert around F4 (MIDI 65) ---
    axis_f4 = note_to_midi('F4') # 65
    notes = [2 * axis_f4 - n for n in notes]
    # Result: [69, 76, 61, 66, 59, 57, 61, 56]
    
    # --- Step 8: Transposed up one Octave (+12 semitones) ---
    notes = [n + 12 for n in notes]
    # Final MIDI Result: [81, 88, 73, 78, 71, 69, 73, 68]
    
    # --- Convert final MIDI list back to note strings and print ---
    final_note_names = [midi_to_note(n) for n in notes]
    
    print(" ".join(final_note_names))

if __name__ == '__main__':
    main()