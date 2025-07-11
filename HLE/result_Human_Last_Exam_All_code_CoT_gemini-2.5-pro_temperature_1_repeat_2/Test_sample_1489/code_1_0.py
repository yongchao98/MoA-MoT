import math

def main():
    """
    Performs a series of musical transformations on a sequence of notes.
    """

    # Helper dictionaries for note-to-MIDI conversion
    NOTE_TO_PC = {'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3, 'E': 4,
                  'F': 5, 'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9,
                  'A#': 10, 'Bb': 10, 'B': 11}
    PC_TO_NOTE = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#',
                  7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}

    def note_to_midi(note_str):
        """Converts a note string (e.g., 'A4') to a MIDI number."""
        name = note_str[:-1]
        octave = int(note_str[-1])
        pc = NOTE_TO_PC[name]
        return 12 * (octave + 1) + pc

    def midi_to_note(midi_num):
        """Converts a MIDI number to a note string (e.g., 'A4')."""
        octave = (midi_num // 12) - 1
        pc = midi_num % 12
        name = PC_TO_NOTE[pc]
        return f"{name}{octave}"

    # Initial sequence
    initial_notes = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # Operation 1: Transpose up Major Second (+2)
    midi_sequence = [n + 2 for n in midi_sequence]

    # Operation 2: Invert Around E4
    pivot_e4 = note_to_midi("E4")
    midi_sequence = [(2 * pivot_e4) - n for n in midi_sequence]

    # Operation 3: Retrograde
    midi_sequence.reverse()

    # Operation 4: Augment intervals by adding 3 semitones
    augmented_sequence = [midi_sequence[0]]
    for i in range(1, len(midi_sequence)):
        original_interval = midi_sequence[i] - midi_sequence[i-1]
        new_note = augmented_sequence[-1] + original_interval + 3
        augmented_sequence.append(new_note)
    midi_sequence = augmented_sequence

    # Operation 5: Change to Dorian Mode Starting from D4
    # D Dorian scale consists of pitch classes {0, 2, 4, 5, 7, 9, 11}
    dorian_pcs = {0, 2, 4, 5, 7, 9, 11}
    quantized_sequence = []
    for n in midi_sequence:
        pc = n % 12
        if pc in dorian_pcs:
            quantized_sequence.append(n)
        else:
            # Find nearest note in scale, rounding up on a tie
            dist_down = 1
            while (pc - dist_down) % 12 not in dorian_pcs:
                dist_down += 1
            dist_up = 1
            while (pc + dist_up) % 12 not in dorian_pcs:
                dist_up += 1
            if dist_up < dist_down:
                quantized_sequence.append(n + dist_up)
            else: # Rounds up on tie
                quantized_sequence.append(n + dist_up)
    midi_sequence = quantized_sequence

    # Operation 6: Transpose down Minor Third (-3)
    midi_sequence = [n - 3 for n in midi_sequence]

    # Operation 7: Invert around F4
    pivot_f4 = note_to_midi("F4")
    midi_sequence = [(2 * pivot_f4) - n for n in midi_sequence]

    # Operation 8: Transpose up one Octave (+12)
    midi_sequence = [n + 12 for n in midi_sequence]

    # Convert final MIDI sequence back to note names
    final_notes = [midi_to_note(n) for n in midi_sequence]
    
    # Print the final result
    print(" ".join(final_notes))

if __name__ == "__main__":
    main()