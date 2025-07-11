import math

def main():
    """
    Performs a series of eight musical operations on a sequence of notes.
    """

    # --- Helper Functions ---
    
    NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    NOTE_MAP = {name: i for i, name in enumerate(NOTE_NAMES)}
    for i, name in enumerate(NOTE_NAMES):
        if len(name) == 2:
            flat_name = NOTE_NAMES[(i-1)%12] + 'b'
            if flat_name.endswith('##b'): # e.g. C##b -> C#
                flat_name = flat_name[:-3] + '#'
            NOTE_MAP[flat_name] = i

    def note_to_midi(note_str):
        """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
        note_str = note_str.strip().upper()
        
        name = note_str[:-1]
        octave_char = note_str[-1]

        if not octave_char.isdigit():
             raise ValueError(f"Invalid note format: {note_str}")
        
        if name in NOTE_MAP:
             pitch_class = NOTE_MAP[name]
        else:
             pitch_class = NOTE_MAP[name[0]]
             if '#' in name:
                pitch_class = (pitch_class + name.count('#')) % 12
             if 'B' in name.upper() and '#' not in name.upper() and 'b' in name:
                pitch_class = (pitch_class - name.count('b')) % 12
        
        octave = int(octave_char)
        return (octave + 1) * 12 + pitch_class

    def midi_to_note(midi_num):
        """Converts a MIDI number to a note string."""
        octave = midi_num // 12 - 1
        pitch_class_index = midi_num % 12
        note_name = NOTE_NAMES[pitch_class_index]
        return f"{note_name}{octave}"

    # --- Initial State ---
    
    # Sequence of eight musical notes
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    current_sequence = [note_to_midi(n) for n in initial_notes_str.split()]
    
    # --- Operation 1: Transpose up Major Second ---
    current_sequence = [n + 2 for n in current_sequence]
    
    # --- Operation 2: Invert Around E4 ---
    axis_e4 = note_to_midi('E4')
    current_sequence = [(2 * axis_e4) - n for n in current_sequence]
    
    # --- Operation 3: Retrograde the sequence ---
    current_sequence = list(reversed(current_sequence))
    
    # --- Operation 4: Augment the intervals ---
    # Formula: new_note[i] = old_note[i] + 3*i
    current_sequence = [n + (3 * i) for i, n in enumerate(current_sequence)]
    
    # --- Operation 5: Change to Dorian Mode Starting from D4 ---
    # D Dorian scale degrees (from C): D(2), E(4), F(5), G(7), A(9), B(11), C(0)
    dorian_degrees = {0, 2, 4, 5, 7, 9, 11}
    d_dorian_notes = [i for i in range(128) if (i % 12) in dorian_degrees]
    
    new_sequence = []
    for note in current_sequence:
        min_dist = min(abs(note - d_note) for d_note in d_dorian_notes)
        # Find all notes with the minimum distance, and choose the lowest pitch in case of a tie
        tied_notes = [d_note for d_note in d_dorian_notes if abs(note - d_note) == min_dist]
        new_sequence.append(min(tied_notes))
    current_sequence = new_sequence

    # --- Operation 6: Transpose down Minor Third ---
    current_sequence = [n - 3 for n in current_sequence]
    
    # --- Operation 7: Invert around F4 ---
    axis_f4 = note_to_midi('F4')
    current_sequence = [(2 * axis_f4) - n for n in current_sequence]

    # --- Operation 8: Transposed up one Octave ---
    current_sequence = [n + 12 for n in current_sequence]

    # --- Final Output ---
    final_note_names = [midi_to_note(n) for n in current_sequence]
    final_answer_str = " ".join(final_note_names)
    
    print(final_answer_str)

if __name__ == "__main__":
    main()