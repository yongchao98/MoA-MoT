def solve_music_sequence():
    """
    Solves the musical note transformation problem.
    """

    # --- Step 1: Helper functions for MIDI conversion ---
    NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    def note_to_midi(note_str):
        """Converts a note string like 'A4' to a MIDI number."""
        note_str = note_str.upper()
        # Handle sharp notes (e.g., 'C#4') vs natural notes (e.g., 'A4')
        if '#' in note_str:
            pitch_class = note_str[:2]
            octave = int(note_str[2:])
        else:
            pitch_class = note_str[0]
            octave = int(note_str[1:])
        
        note_index = NOTE_NAMES.index(pitch_class)
        # MIDI standard: C4 = 60, A4 = 69
        # Formula: midi = 12 * (octave + 1) + note_index
        return 12 * (octave + 1) + note_index

    def midi_to_note(midi_num):
        """Converts a MIDI number back to a note string."""
        octave = (midi_num // 12) - 1
        note_index = midi_num % 12
        pitch_class = NOTE_NAMES[note_index]
        return f"{pitch_class}{octave}"

    # --- Step 2: Initial sequence ---
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    notes = [note_to_midi(n) for n in initial_notes_str]

    # --- Step 3: Apply the 8 operations sequentially ---

    # 1. Transpose up Major Second (+2 semitones)
    notes = [n + 2 for n in notes]

    # 2. Invert Around E4
    axis_e4 = note_to_midi('E4')
    notes = [(2 * axis_e4) - n for n in notes]

    # 3. Retrograde the sequence
    notes.reverse()

    # 4. Augment intervals by adding 3 semitones
    if len(notes) > 1:
        augmented_notes = [notes[0]]
        for i in range(1, len(notes)):
            original_interval = notes[i] - notes[i-1]
            new_interval = original_interval + 3
            next_note = augmented_notes[-1] + new_interval
            augmented_notes.append(next_note)
        notes = augmented_notes

    # 5. Change to Dorian Mode Starting from D4
    dorian_root_midi = note_to_midi('D4')
    # Dorian scale steps from root: [0, 2, 3, 5, 7, 9, 10]
    dorian_steps = {0, 2, 3, 5, 7, 9, 10}
    dorian_scale_notes = [n for n in range(128) if (n - dorian_root_midi) % 12 in dorian_steps]
    
    mapped_notes = []
    for note in notes:
        # Find the closest note in the Dorian scale.
        # The key sorts first by distance, then by value to handle ties (choosing the lower note).
        closest_note = min(dorian_scale_notes, key=lambda d_note: (abs(note - d_note), d_note))
        mapped_notes.append(closest_note)
    notes = mapped_notes

    # 6. Transpose down Minor Third (-3 semitones)
    notes = [n - 3 for n in notes]

    # 7. Invert around F4
    axis_f4 = note_to_midi('F4')
    notes = [(2 * axis_f4) - n for n in notes]

    # 8. Transposed up one Octave (+12 semitones)
    notes = [n + 12 for n in notes]

    # --- Step 4 & 5: Final Conversion and Output ---
    final_note_names = [midi_to_note(n) for n in notes]
    print(" ".join(final_note_names))

solve_music_sequence()