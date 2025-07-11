import math

def solve_music_sequence():
    """
    Solves the entire music transformation problem.
    """

    # --- Helper Functions ---
    NOTE_PITCH_CLASSES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]

    def note_to_midi(note_str: str) -> int:
        """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
        note_name = note_str[:-1]
        octave = int(note_str[-1])
        pitch_class = NOTE_PITCH_CLASSES.index(note_name)
        return (octave + 1) * 12 + pitch_class

    def midi_to_note(midi_num: int) -> str:
        """Converts a MIDI number to a note string like 'A4'."""
        octave = midi_num // 12 - 1
        note_index = midi_num % 12
        note_name = NOTE_PITCH_CLASSES[note_index]
        return f"{note_name}{octave}"

    # --- Transformation Functions ---

    def transpose(notes, semitones):
        return [note + semitones for note in notes]

    def invert(notes, axis_note_str):
        axis_midi = note_to_midi(axis_note_str)
        return [2 * axis_midi - note for note in notes]

    def retrograde(notes):
        return notes[::-1]

    def augment_intervals(notes, semitones_to_add):
        if not notes:
            return []
        new_notes = [notes[0]]
        for i in range(len(notes) - 1):
            interval = notes[i+1] - notes[i]
            augmented_interval = interval + semitones_to_add
            new_notes.append(new_notes[-1] + augmented_interval)
        return new_notes

    def map_to_dorian(notes, root_note_str):
        root_midi = note_to_midi(root_note_str)
        root_pc = root_midi % 12
        dorian_intervals = [0, 2, 3, 5, 7, 9, 10] # W-H-W-W-W-H-W
        dorian_pcs = sorted([(root_pc + i) % 12 for i in dorian_intervals])
        
        mapped_notes = []
        for note in notes:
            pc = note % 12
            
            # Find the closest pitch class in the Dorian scale
            min_dist = float('inf')
            best_pcs = []
            for scale_pc in dorian_pcs:
                dist = abs(pc - scale_pc)
                dist = min(dist, 12 - dist) # Account for octave wrapping
                if dist < min_dist:
                    min_dist = dist
                    best_pcs = [scale_pc]
                elif dist == min_dist:
                    best_pcs.append(scale_pc)
            
            # Tie-breaking rule: round down (choose the lower pitch class)
            best_pc = min(best_pcs)

            # Reconstruct the note in the closest octave
            octave_base = note - pc
            candidates = [octave_base + best_pc, octave_base - 12 + best_pc, octave_base + 12 + best_pc]
            closest_note = min(candidates, key=lambda x: abs(x - note))
            mapped_notes.append(closest_note)
            
        return mapped_notes

    # --- Main Execution ---

    # Initial sequence: A4 C5 F4 E4 G4 C4 B4 D4
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    current_sequence = [note_to_midi(n) for n in initial_notes_str]

    # 1. Transpose up Major Second (+2 semitones)
    current_sequence = transpose(current_sequence, 2)

    # 2. Invert Around E4
    current_sequence = invert(current_sequence, 'E4')

    # 3. Retrograde the sequence
    current_sequence = retrograde(current_sequence)

    # 4. Augment intervals by adding 3 semitones
    current_sequence = augment_intervals(current_sequence, 3)

    # 5. Change to Dorian Mode Starting from D4
    current_sequence = map_to_dorian(current_sequence, 'D4')

    # 6. Transpose down Minor Third (-3 semitones)
    current_sequence = transpose(current_sequence, -3)

    # 7. Invert around F4
    current_sequence = invert(current_sequence, 'F4')

    # 8. Transposed up one Octave (+12 semitones)
    current_sequence = transpose(current_sequence, 12)

    # Convert final MIDI numbers back to note names
    final_note_names = [midi_to_note(n) for n in current_sequence]

    # Print the final result
    print(" ".join(final_note_names))

solve_music_sequence()