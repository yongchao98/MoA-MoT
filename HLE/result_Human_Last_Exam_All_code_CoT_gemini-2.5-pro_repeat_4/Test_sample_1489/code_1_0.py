import re
from bisect import bisect_left

def solve_music_sequence():
    """
    Solves the entire music theory problem by applying a sequence of operations to a list of notes.
    """
    # --- Step 1: Define Note Representation Helpers ---

    NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    NOTE_MAP = {name: i for i, name in enumerate(NOTE_NAMES)}

    def note_to_midi(note_str: str) -> int:
        """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
        match = re.match(r"([A-G])(#?b?)(\d+)", note_str)
        if not match:
            raise ValueError(f"Invalid note format: {note_str}")
        
        name, accidental, octave_str = match.groups()
        note_index = NOTE_MAP[name + accidental]
        octave = int(octave_str)
        
        # MIDI standard: C4 is 60. Formula: 12 * (octave + 1) + note_index
        return 12 * (octave + 1) + note_index

    def midi_to_note(midi_num: int) -> str:
        """Converts a MIDI number to a note string."""
        if not (0 <= midi_num <= 127):
            raise ValueError("MIDI number out of range (0-127)")
        
        octave = (midi_num // 12) - 1
        note_index = midi_num % 12
        note_name = NOTE_NAMES[note_index]
        return f"{note_name}{octave}"

    # --- Step 2: Define Transformation Functions ---

    def transpose(notes, semitones):
        return [n + semitones for n in notes]

    def invert(notes, axis_note_str):
        axis_midi = note_to_midi(axis_note_str)
        return [2 * axis_midi - n for n in notes]

    def retrograde(notes):
        return notes[::-1]

    def augment_intervals(notes, semitones_to_add):
        if not notes:
            return []
        # The k-th note is augmented by (k-1) * semitones_to_add
        return [notes[i] + i * semitones_to_add for i in range(len(notes))]

    def map_to_dorian_mode(notes, root_note_str):
        # D Dorian pitch classes (relative to C=0): D, E, F, G, A, B, C
        # MIDI pitch classes: 2, 4, 5, 7, 9, 11, 0
        scale_pcs = sorted([0, 2, 4, 5, 7, 9, 11])
        
        def find_closest_pc(target_pc):
            # Find the insertion point for the target in the sorted list of scale tones
            pos = bisect_left(scale_pcs, target_pc)
            if pos == 0:
                # Closer to scale_pcs[0] or scale_pcs[-1] - 12?
                if abs(target_pc - scale_pcs[0]) <= abs(target_pc - (scale_pcs[-1] - 12)):
                    return scale_pcs[0]
                else:
                    return scale_pcs[-1]
            if pos == len(scale_pcs):
                 # Closer to scale_pcs[-1] or scale_pcs[0] + 12?
                if abs(target_pc - scale_pcs[-1]) <= abs(target_pc - (scale_pcs[0] + 12)):
                    return scale_pcs[-1]
                else:
                    return scale_pcs[0]

            before = scale_pcs[pos - 1]
            after = scale_pcs[pos]
            
            # Tie-break rule: round down
            if after - target_pc < target_pc - before:
                return after
            else:
                return before

        mapped_notes = []
        for midi_note in notes:
            original_pc = midi_note % 12
            if original_pc not in scale_pcs:
                closest_pc = find_closest_pc(original_pc)
                # Adjust original note by the difference, preserving octave context
                diff = closest_pc - original_pc
                # Handle wrapping around the octave (e.g. C -> B is -1, not +11)
                if diff > 6: diff -= 12
                if diff < -6: diff += 12
                mapped_notes.append(midi_note + diff)
            else:
                mapped_notes.append(midi_note)
        return mapped_notes

    # --- Step 3: Execute the sequence of operations ---

    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    current_notes = [note_to_midi(n) for n in initial_notes_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    current_notes = transpose(current_notes, 2)
    
    # 2. Invert Around E4
    current_notes = invert(current_notes, "E4")
    
    # 3. Retrograde the sequence
    current_notes = retrograde(current_notes)
    
    # 4. Augment intervals by adding 3 semitones
    current_notes = augment_intervals(current_notes, 3)
    
    # 5. Change to Dorian Mode Starting from D4
    current_notes = map_to_dorian_mode(current_notes, "D4")
    
    # 6. Transpose down Minor Third (-3 semitones)
    current_notes = transpose(current_notes, -3)
    
    # 7. Invert around F4
    current_notes = invert(current_notes, "F4")
    
    # 8. Transposed up one Octave (+12 semitones)
    current_notes = transpose(current_notes, 12)

    # --- Step 4: Convert final MIDI list to string and print ---
    final_note_names = [midi_to_note(n) for n in current_notes]
    print(" ".join(final_note_names))

solve_music_sequence()