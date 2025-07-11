import math

def solve_music_sequence():
    """
    Solves the music theory problem by applying a series of transformations to a sequence of notes.
    """

    # --- Step 1: Setup and Helper Functions ---

    NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    NOTE_MAP = {name: i for i, name in enumerate(NOTE_NAMES)}

    def note_to_midi(note_str):
        """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
        name = note_str[:-1]
        octave = int(note_str[-1])
        if len(name) > 1 and name[1] in '#b':
            octave = int(note_str[2:])
            name = note_str[:2]
        else:
            octave = int(note_str[1:])
            name = note_str[:1]
        
        # Handle flat names for robustness, though not in the input
        flat_map = {"Db": "C#", "Eb": "D#", "Gb": "F#", "Ab": "G#", "Bb": "A#"}
        name = flat_map.get(name, name)

        return NOTE_MAP[name] + (octave + 1) * 12

    def midi_to_note(midi_num):
        """Converts a MIDI number back to a note string like 'A4'."""
        if not isinstance(midi_num, int):
            return "Invalid"
        octave = midi_num // 12 - 1
        note_index = midi_num % 12
        return f"{NOTE_NAMES[note_index]}{octave}"

    # --- Step 2: Define Transformation Functions ---

    def transpose(notes, semitones):
        return [note + semitones for note in notes]

    def invert(notes, axis_midi):
        return [2 * axis_midi - note for note in notes]

    def retrograde(notes):
        return notes[::-1]

    def augment_intervals(notes, amount):
        # new_note[i] = old_note[i] + amount * i
        return [notes[i] + amount * i for i in range(len(notes))]

    def snap_to_mode(notes, root_midi, mode_intervals):
        """Maps each note to the closest note in the specified mode."""
        scale_notes = []
        # Generate scale notes over a wide range of octaves
        for i in range(-5, 6):
            octave_root = root_midi + 12 * i
            for interval in mode_intervals:
                scale_notes.append(octave_root + interval)
        scale_notes = sorted(list(set(scale_notes)))

        new_notes = []
        for note in notes:
            # Find the closest note, rounding up on a tie
            min_dist = float('inf')
            closest_note = None
            for scale_note in scale_notes:
                dist = abs(note - scale_note)
                if dist < min_dist:
                    min_dist = dist
                    closest_note = scale_note
                elif dist == min_dist:
                    closest_note = max(closest_note, scale_note)
            new_notes.append(closest_note)
        return new_notes

    # --- Step 3: Apply Operations Sequentially ---

    # Initial sequence
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    notes = [note_to_midi(n) for n in initial_notes_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    notes = transpose(notes, 2)

    # 2. Invert Around E4 (MIDI 64)
    notes = invert(notes, note_to_midi("E4"))

    # 3. Retrograde the sequence
    notes = retrograde(notes)

    # 4. Augment intervals by adding 3 semitones
    notes = augment_intervals(notes, 3)

    # 5. Change to Dorian Mode Starting from D4 (root 62)
    d4_midi = note_to_midi("D4")
    dorian_intervals = [0, 2, 3, 5, 7, 9, 10] # T-S-T-T-T-S-T
    notes = snap_to_mode(notes, d4_midi, dorian_intervals)

    # 6. Transpose down Minor Third (-3 semitones)
    notes = transpose(notes, -3)

    # 7. Invert around F4 (MIDI 65)
    notes = invert(notes, note_to_midi("F4"))

    # 8. Transposed up one Octave (+12 semitones)
    notes = transpose(notes, 12)

    # --- Step 4: Final Conversion and Output ---
    final_note_names = [midi_to_note(n) for n in notes]
    print(" ".join(final_note_names))

solve_music_sequence()