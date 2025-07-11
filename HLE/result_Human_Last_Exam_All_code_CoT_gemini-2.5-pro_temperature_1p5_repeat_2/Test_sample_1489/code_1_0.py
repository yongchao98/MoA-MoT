import math

def solve_music_sequence():
    """
    Applies a series of eight musical operations to a sequence of notes
    and prints the final result.
    """

    # Helper functions for converting between note names and MIDI numbers
    NOTE_MAP = {'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3, 'E': 4, 'F': 5,
                'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11}
    NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    def note_to_midi(note_str):
        """Converts a note string like 'A4' to a MIDI number."""
        # Handles note names like 'C#' or 'Db'
        if len(note_str) > 2 and note_str[1] in '#b':
            name = note_str[:2]
            octave = int(note_str[2:])
        else:
            name = note_str[:1]
            octave = int(note_str[1:])
        return NOTE_MAP[name] + (octave + 1) * 12

    def midi_to_note(midi_num):
        """Converts a MIDI number to a note string like 'A4'."""
        octave = midi_num // 12 - 1
        note_index = midi_num % 12
        return f"{NOTE_NAMES[note_index]}{octave}"

    # ------------------
    # Step 0: Initial Sequence
    # ------------------
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    notes = [note_to_midi(n) for n in initial_notes_str.split()]
    # Initial MIDI: [69, 72, 65, 64, 67, 60, 71, 62]

    # ------------------
    # Operation 1: Transpose up Major Second (+2 semitones)
    # ------------------
    notes = [n + 2 for n in notes]
    # Result: [71, 74, 67, 66, 69, 62, 73, 64]

    # ------------------
    # Operation 2: Invert Around E4 (MIDI 64)
    # ------------------
    axis_e4 = note_to_midi('E4')
    notes = [(2 * axis_e4) - n for n in notes]
    # Result: [57, 54, 61, 62, 59, 66, 55, 64]

    # ------------------
    # Operation 3: Retrograde the sequence
    # ------------------
    notes.reverse()
    # Result: [64, 55, 66, 59, 62, 61, 54, 57]

    # ------------------
    # Operation 4: Augment the intervals between consecutive notes by adding 3 semitones.
    # ------------------
    augmented_notes = []
    if notes:
        augmented_notes.append(notes[0])
        for i in range(1, len(notes)):
            original_interval = notes[i] - notes[i-1]
            new_interval = original_interval + 3
            new_note = augmented_notes[i-1] + new_interval
            augmented_notes.append(new_note)
        notes = augmented_notes
    # Result: [64, 58, 72, 68, 74, 76, 72, 78]

    # ------------------
    # Operation 5: Change to Dorian Mode Starting from D4 (Quantize to D Dorian scale)
    # ------------------
    d_dorian_intervals = [0, 2, 3, 5, 7, 9, 10]
    d_dorian_notes = []
    # Generate scale notes over a wide range of octaves to find the closest match
    for octave_offset in range(-1, 10):
        root = 12 * (octave_offset + 1) + NOTE_MAP['D']
        for interval in d_dorian_intervals:
            d_dorian_notes.append(root + interval)
    d_dorian_notes.sort()

    quantized_notes = []
    for note in notes:
        # Find the closest note in the scale.
        # On a tie in distance, the one with the smaller MIDI number (lower pitch) is chosen.
        closest_note = min(d_dorian_notes, key=lambda x: (abs(x - note), x))
        quantized_notes.append(closest_note)
    notes = quantized_notes
    # Result: [64, 57, 72, 67, 74, 76, 72, 77]

    # ------------------
    # Operation 6: Transpose down Minor Third (-3 semitones)
    # ------------------
    notes = [n - 3 for n in notes]
    # Result: [61, 54, 69, 64, 71, 73, 69, 74]

    # ------------------
    # Operation 7: Invert around F4 (MIDI 65)
    # ------------------
    axis_f4 = note_to_midi('F4')
    notes = [(2 * axis_f4) - n for n in notes]
    # Result: [69, 76, 61, 66, 59, 57, 61, 56]

    # ------------------
    # Operation 8: Transposed up one Octave (+12 semitones)
    # ------------------
    notes = [n + 12 for n in notes]
    # Final MIDI: [81, 88, 73, 78, 71, 69, 73, 68]

    # ------------------
    # Final conversion and output
    # ------------------
    final_note_names = [midi_to_note(n) for n in notes]
    print(' '.join(final_note_names))

solve_music_sequence()